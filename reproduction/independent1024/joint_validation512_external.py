"""Fixed-output external/structural sensitivity for all seven raw512 policies.

prepare records coverage and missing fold predictions without loading models.
run uses existing published CombiMOTS checkpoints on CPU, missing folds only.
analyze requires complete cached predictions and never loads a property model.
self-test uses injected chemistry/predictions and temporary mock files only.

No task RF, Graph search, candidate selection, training, or deployment occurs
here. Native/Graph execution unknowns are bounded rather than set to failure.
One writer per output directory; all fold batches are saved for restart.
"""
import argparse
from collections import Counter
import csv
import importlib.metadata
import json
import math
import os
from pathlib import Path
import statistics
import time

from joint_validation512_analysis import CONTRASTS, exact_stratified_bootstrap, limits, percentile
from joint_validation512_graph import POLICIES, THRESHOLDS, identity, validate_sources

DEV = Path('/home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development')
TASKS = ('gsk3b_activity', 'jnk3_activity')
TASK_INDEX = {'gsk3b': 0, 'jnk3': 1}
EXTERNAL_THRESHOLDS = (.3, .5, .7)
MODEL_REFERENCE = {
    'repository': 'https://github.com/Tibogoss/CombiMOTS',
    'published_commit': '3f644b727982a95b4f87115ccb8d6b4c918a02ba',
    'software': 'Chemprop1.5.2', 'tasks': list(TASKS), 'models': 10,
}
FILTERS = ('all', 'no_pains', 'no_pains_brenk', 'mw_le500', 'source_qed_tolerance')
LIMITATIONS = [
    'Same broad ExCAPE data source as optimized RF; cross-architecture sensitivity, '
    'not independent biological validation or measured potency.',
    'Checkpoint args reference data/excape+aixfuse_curated.csv, train_data_size77581; '
    'public CSV95159rows does not identify complete original training rows/order.',
    'No reconstructed OOF, complete out-of-training or scaffold exclusion claim. '
    'Public-table exact matches are descriptive provenance only.',
    'Probability thresholds are uncalibrated outputs. Fold SD is model dispersion, '
    'not biological uncertainty.',
    'Only eight training runs, four within each of two fixed task predictors. '
    'Conditional run-level bootstrap; individual slots are not independent replicates.',
    'Known fallback structure descriptors/probabilities remain operational observations. '
    'Execution-unknown policy replacements may differ; bounds are not confidence intervals.',
]


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def valid_vector(value):
    return (isinstance(value, list) and len(value) == 2 and all(
        isinstance(v, (int, float)) and not isinstance(v, bool) and
        math.isfinite(v) and 0 <= v <= 1 for v in value))


def configure_cpu():
    os.environ['CUDA_VISIBLE_DEVICES'] = ''
    for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[variable] = '2'


class Chemistry:
    """Lazy RDKit descriptors, matching old256 radius2/2048/no-chirality settings."""
    def __init__(self):
        from rdkit import Chem, DataStructs, RDLogger, rdBase
        from rdkit.Chem import AllChem, Descriptors, QED
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
        from rdkit.Chem.Scaffolds import MurckoScaffold
        RDLogger.DisableLog('rdApp.*')
        self.Chem, self.DataStructs = Chem, DataStructs
        self.AllChem, self.Descriptors, self.QED = AllChem, Descriptors, QED
        self.MurckoScaffold = MurckoScaffold
        self.version = rdBase.rdkitVersion
        self.molecules, self.fingerprints = {}, {}
        self.catalogs = {}
        for name, parts in [('pains', ('PAINS_A', 'PAINS_B', 'PAINS_C')), ('brenk', ('BRENK',))]:
            parameters = FilterCatalogParams()
            for part in parts:
                parameters.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, part))
            self.catalogs[name] = FilterCatalog(parameters)

    def molecule(self, smiles):
        if not smiles:
            return None
        if smiles not in self.molecules:
            self.molecules[smiles] = self.Chem.MolFromSmiles(smiles)
        return self.molecules[smiles]

    def canonical(self, smiles):
        molecule = self.molecule(smiles)
        return (self.Chem.MolToSmiles(molecule, isomericSmiles=True)
                if molecule is not None else None)

    def profile(self, smiles):
        molecule = self.molecule(smiles)
        if molecule is None:
            raise ValueError('A recorded valid output cannot be parsed')
        return {
            'canonical_smiles': self.canonical(smiles), 'qed': float(self.QED.qed(molecule)),
            'molecular_weight': float(self.Descriptors.MolWt(molecule)),
            'heavy_atoms': int(molecule.GetNumHeavyAtoms()),
            'components': len(self.Chem.GetMolFrags(molecule)),
            'murcko_scaffold': self.MurckoScaffold.MurckoScaffoldSmiles(
                mol=molecule, includeChirality=False),
            **{name: sorted({entry.GetDescription() for entry in catalog.GetMatches(molecule)})
               for name, catalog in self.catalogs.items()},
        }

    def similarity(self, left, right):
        if not left or not right:
            return None
        for smiles in (left, right):
            if smiles not in self.fingerprints:
                self.fingerprints[smiles] = self.AllChem.GetMorganFingerprintAsBitVect(
                    self.molecule(smiles), 2, nBits=2048)
        return float(self.DataStructs.TanimotoSimilarity(
            self.fingerprints[left], self.fingerprints[right]))


def fixed_inputs(study, reports, chemistry):
    design = read(study / 'design.json')
    sources = read(study / 'sources.json')['rows']
    validate_sources(sources, design)
    if len(design['runs']) != 8 or design['raw_denominator'] != 512:
        raise ValueError('Expected the recorded eight-run raw512 study')
    planned = [p['name'] for p in design['policies']]
    if set(planned) != set(POLICIES) or len(planned) != len(POLICIES):
        raise ValueError('Expected all seven prespecified policies')
    source_map = {s['sample_id']: s for s in sources}
    observations = [r for path in reports for r in read(path)['rows']]
    expected = {(s['sample_id'], p) for s in sources for p in POLICIES}
    observed = [(r['sample_id'], r['policy']) for r in observations]
    if set(observed) != expected or len(observed) != len(expected):
        raise ValueError('All512 slots and all7 policies must occur exactly once')
    by_key = {(r['sample_id'], r['policy']): r for r in observations}
    rows = [by_key[(s['sample_id'], p)] for s in sources for p in POLICIES]
    for source in sources:
        smiles = source.get('canonical_smiles')
        if source['valid'] != bool(smiles):
            raise ValueError('Source validity/canonical structure disagree')
        if smiles and chemistry.canonical(smiles) != smiles:
            raise ValueError('Source canonical identity differs; do not alter recorded sources')
    for row in rows:
        source = source_map[row['sample_id']]
        if row['source_identity'] != identity(source):
            raise ValueError('Graph policy source identity differs')
        for key, source_key in [('source_smiles', 'canonical_smiles'),
                                ('source_score', 'task_score'), ('source_qed', 'qed')]:
            if row[key] != source.get(source_key):
                raise ValueError('Policy source field differs: ' + key)
        for key in ('task', 'replicate', 'cohort_id', 'position'):
            if row[key] != source[key]:
                raise ValueError('Policy slot identity differs: ' + key)
        if row['status'] in ('processing', 'pending_candidate_stage'):
            raise ValueError('Candidate policies are still pending; no fixed-output analysis yet')
        for endpoint in ('solved', 'joint'):
            limits(row, endpoint)
        final = row['final_smiles']
        if final and chemistry.canonical(final) != final:
            raise ValueError('Final output is not recorded canonical SMILES; do not reselect')
        if not final and (row['solved'] is True or row['joint'] is True):
            raise ValueError('Missing final structure cannot be scientifically successful')
        if not source['valid'] and (final or not row['scientific_endpoint_available']):
            raise ValueError('Invalid source must retain its known invalid output')
        if not row['scientific_endpoint_available'] and final != source['canonical_smiles']:
            raise ValueError('Unknown policy must retain the recorded source fallback')
    return design, sources, rows


def default_history(dev):
    paths = [
        dev / 'external_property_raw96_v1/predictions.json',
        dev / 'external_property_raw96_v1/predictions_with_depth3.json',
        dev / 'external_property_independent256_v1/predictions.json',
    ]
    for fold in range(10):
        paths.extend([
            dev / f'external_property_raw96_v1/fold_{fold:02d}_predictions.json',
            dev / f'external_property_raw96_v1/depth3_fold_{fold:02d}_new_predictions.json',
            dev / f'external_property_independent256_v1/new_fold_{fold:02d}_predictions.json',
            dev / f'budget_replay256_v1/external_new_fold_{fold:02d}.json',
        ])
    return sorted({p for p in paths if p.exists()}, key=str)


def cache_fold(path, data):
    value = data.get('fold')
    if not isinstance(value, int) or isinstance(value, bool) or not 0 <= value < 10:
        raise ValueError('Expected fold index0..9: ' + str(path))
    return value


def cache_pairs(path, data):
    predictions = data.get('predictions')
    if isinstance(predictions, dict):
        return predictions.items()
    if isinstance(predictions, list):
        smiles = data.get('smiles', [])
        if len(smiles) != len(predictions) or len(set(smiles)) != len(smiles):
            raise ValueError('Malformed legacy parallel smiles/prediction arrays')
        return zip(smiles, predictions)
    raise ValueError('Missing predictions: ' + str(path))


def load_history(paths, canonical):
    """Merge every compatible cache, retaining fold-level provenance.

    File order is lexical, independent of scores/outcomes. Legacy no-task-label
    fold files are accepted only under their documented old experiment paths.
    """
    values = [{} for _ in range(10)]
    provenance = [{} for _ in range(10)]
    files, duplicate_entries = [], 0
    for path in sorted(set(map(Path, paths)), key=str):
        data = read(path)
        if data.get('tasks') is not None and data['tasks'] != list(TASKS):
            raise ValueError('Cache task-column order differs: ' + str(path))
        if data.get('models') is not None and data['models'] != 10:
            raise ValueError('Cache model count differs')
        reference = data.get('model_reference')
        if reference is not None and reference != MODEL_REFERENCE:
            raise ValueError('Cache records a different property model')
        if data.get('tasks') is None:
            documented = path.parent.name in ('external_property_raw96_v1', 'budget_replay256_v1')
            if not documented or data.get('fold') is None:
                raise ValueError('Unlabelled cache is not a documented legacy fold file')
        added = reused = 0
        for raw, record in cache_pairs(path, data):
            smiles = canonical(raw)
            if not smiles:
                raise ValueError('Invalid structure in historical property cache')
            if isinstance(record, dict):
                folds = record.get('folds')
                if not isinstance(folds, list) or len(folds) != 10 or not all(map(valid_vector, folds)):
                    raise ValueError('Historical mean is insufficient: require all10 fold vectors')
                mean = record.get('mean')
                if not valid_vector(mean) or any(abs(mean[t] - statistics.mean(
                        v[t] for v in folds)) > 1e-7 for t in range(2)):
                    raise ValueError('Historical mean/folds disagree')
                requested = enumerate(folds)
            else:
                requested = [(cache_fold(path, data), record)]
            for fold, vector in requested:
                if not valid_vector(vector):
                    raise ValueError('Malformed property vector: ' + str(path))
                previous = values[fold].get(smiles)
                if previous is not None:
                    if any(abs(a - b) > 1e-7 for a, b in zip(previous, vector)):
                        raise ValueError('Conflicting same-model/canonical/fold cache predictions')
                    reused += 1
                    duplicate_entries += 1
                    provenance[fold][smiles]['equivalent_cache_paths'].append(str(path))
                else:
                    values[fold][smiles] = list(vector)
                    provenance[fold][smiles] = {
                        'selected_path': str(path), 'recorded_smiles': raw,
                        'equivalent_cache_paths': [], 'kind': 'historical',
                    }
                    added += 1
        files.append({'path': str(path), 'new_fold_structure_entries': added,
                      'equivalent_fold_structure_entries': reused})
    return values, provenance, {
        'cache_files': files, 'equivalent_entries': duplicate_entries,
        'historical_unique_structures': len(set().union(*(set(v) for v in values))),
        'historical_complete10_structures': len(set.intersection(*(set(v) for v in values))),
        'model_reference': MODEL_REFERENCE,
    }


class PublishedPredictor:
    def __init__(self, checkpoint):
        configure_cpu()
        import torch
        from chemprop.utils import load_args, load_checkpoint
        if importlib.metadata.version('chemprop') != '1.5.2':
            raise RuntimeError('Use the existing Chemprop1.5.2 environment for released checkpoints')
        args = load_args(str(checkpoint))
        if list(args.task_names) != list(TASKS):
            raise ValueError('Published checkpoint task names/order differ')
        torch.set_num_threads(2)
        self.torch = torch
        self.model = load_checkpoint(str(checkpoint), device=torch.device('cpu')).eval()
        if next(self.model.parameters()).device.type != 'cpu':
            raise RuntimeError('Expected CPU-only property predictor')

    def __call__(self, smiles):
        with self.torch.inference_mode():
            values = self.model(batch=[[s] for s in smiles], features_batch=None).cpu().numpy().tolist()
        if len(values) != len(smiles) or not all(map(valid_vector, values)):
            raise ValueError('Unexpected published checkpoint predictions')
        return values


def required_structures(sources, rows):
    structures = {s['canonical_smiles'] for s in sources if s['valid']}
    structures.update(r['final_smiles'] for r in rows if r['final_smiles'])
    pairs = {(s['task'], s['canonical_smiles']) for s in sources if s['valid']}
    pairs.update((r['task'], r['final_smiles']) for r in rows if r['final_smiles'])
    return sorted(structures), sorted(pairs)


def new_fold(out, fold):
    path = out / f'new_fold_{fold:02d}_predictions.json'
    data = read(path) if path.exists() else {
        'fold': fold, 'tasks': list(TASKS), 'model_reference': MODEL_REFERENCE,
        'predictions': {}, 'batches': [],
    }
    if (data.get('fold') != fold or data.get('tasks') != list(TASKS) or
            data.get('model_reference') != MODEL_REFERENCE or
            any(not valid_vector(v) for v in data['predictions'].values())):
        raise ValueError('Malformed current fold cache')
    return path, data


def obtain_predictions(required, history, provenance, assets, out, allow_inference,
                       predictor_factory=PublishedPredictor):
    stored, new_calls, loaded = [], 0, 0
    resumed_entries = historical_entries = 0
    for fold in range(10):
        path, current = new_fold(out, fold)
        for smiles, vector in current['predictions'].items():
            previous = history[fold].get(smiles)
            if previous is not None and any(abs(a - b) > 1e-7 for a, b in zip(previous, vector)):
                raise ValueError('Current/history prediction conflict')
        missing = [s for s in required if s not in history[fold] and s not in current['predictions']]
        if missing and not allow_inference:
            raise RuntimeError(f'Fold{fold} missing{len(missing)} vectors; analyze never loads models')
        if missing:
            checkpoint = assets / f'published/models/gsk3b_jnk3/fold_{fold}/model_0/model.pt'
            predictor = predictor_factory(checkpoint)
            loaded += 1
            try:
                for start in range(0, len(missing), 32):
                    batch = missing[start:start + 32]
                    begin = time.monotonic()
                    vectors = predictor(batch)
                    if len(vectors) != len(batch) or not all(map(valid_vector, vectors)):
                        raise ValueError('Invalid inference batch')
                    current['predictions'].update(zip(batch, vectors))
                    current['batches'].append({
                        'smiles': batch, 'structures': len(batch),
                        'seconds': time.monotonic() - begin,
                    })
                    save(path, current)
                    new_calls += len(batch)
            finally:
                del predictor
        combined = {}
        for smiles in required:
            if smiles in history[fold]:
                combined[smiles] = (history[fold][smiles], provenance[fold][smiles])
                historical_entries += 1
            else:
                combined[smiles] = (current['predictions'][smiles], {
                    'kind': 'current_missing_only', 'selected_path': str(path),
                    'recorded_smiles': smiles, 'equivalent_cache_paths': [],
                })
                if smiles not in missing:
                    resumed_entries += 1
        stored.append(combined)
    predictions = {}
    for smiles in required:
        vectors = [v[smiles][0] for v in stored]
        predictions[smiles] = {
            'mean': [statistics.mean(v[t] for v in vectors) for t in range(2)],
            'fold_sd': [statistics.pstdev(v[t] for v in vectors) for t in range(2)],
            'folds': vectors, 'fold_provenance': [v[smiles][1] for v in stored],
        }
    audit = {
        'required_unique_structures_including_sources': len(required),
        'logical_checkpoint_structure_predictions': len(required) * 10,
        'reused_historical_fold_structure_predictions': historical_entries,
        'reused_current_fold_structure_predictions': resumed_entries,
        'model_structure_predictions_this_invocation': new_calls,
        'loaded_models_this_invocation': loaded, 'historical_recomputed': 0,
        'historical_complete10_in_required': sum(all(s in h for h in history) for s in required),
    }
    return predictions, audit


def public_index(dev, assets, out, chemistry):
    for path in (out / 'public_canonical_index.json',
                 dev / 'external_property_independent256_v1/public_canonical_index.json',
                 assets / 'public_canonical_index.json'):
        if path.exists():
            return read(path), str(path)
    source = assets / 'published/data/GSK3B_JNK3.csv'
    if not source.exists():
        return None, None
    lookup, count, invalid = {}, 0, 0
    with source.open(encoding='utf-8', newline='') as stream:
        for i, row in enumerate(csv.DictReader(stream)):
            count += 1
            smiles = chemistry.canonical(row['smiles'])
            if smiles:
                lookup.setdefault(smiles, []).append({
                    'public_row': i, **{task: row[task] for task in TASKS},
                })
            else:
                invalid += 1
    result = {'source': str(source), 'lookup': lookup, 'public_rows': count,
              'public_unique_valid': len(lookup), 'public_invalid_rows': invalid}
    path = out / 'public_canonical_index.json'
    save(path, result)
    return result, str(path)


def descriptor_records(required, rows, chemistry, out):
    path = out / 'structure_descriptors.json'
    cached = read(path) if path.exists() else {
        'rdkit_version': chemistry.version, 'fingerprint': 'Morgan radius2/2048,chiral=False',
        'scaffold': 'Murcko includeChirality=False; empty acyclic counted separately',
        'profiles': {}, 'source_final_similarities': [],
    }
    if cached['rdkit_version'] != chemistry.version:
        raise ValueError('Descriptor cache RDKit version differs')
    profiles = cached['profiles']
    computed = 0
    for smiles in required:
        if smiles not in profiles:
            profiles[smiles] = chemistry.profile(smiles)
            computed += 1
            if computed % 32 == 0:
                save(path, cached)
    similarities = {(r['source'], r['final']): r['similarity']
                    for r in cached['source_final_similarities']}
    pairs = {(r['source_smiles'], r['final_smiles']) for r in rows}
    for source, final in sorted(pairs, key=lambda x: (str(x[0]), str(x[1]))):
        if (source, final) not in similarities:
            similarities[(source, final)] = chemistry.similarity(source, final)
    cached['source_final_similarities'] = [
        {'source': a, 'final': b, 'similarity': value}
        for (a, b), value in sorted(similarities.items(), key=lambda x: str(x[0]))]
    save(path, cached)
    return profiles, similarities, {'new_structure_profiles': computed,
                                   'required_unique_structures': len(required)}


def enrich(rows, predictions, profiles, similarities, public):
    lookup = public['lookup'] if public is not None else {}
    enriched = []
    for original in rows:
        row = dict(original)
        source, final, index = row['source_smiles'], row['final_smiles'], TASK_INDEX[row['task']]
        probability = predictions[final]['mean'][index] if final else 0.
        source_probability = predictions[source]['mean'][index] if source else 0.
        matches = lookup.get(final, []) if final else []
        route_operational = (row['solved'] if row['scientific_endpoint_available']
                             else row['fallback_solved'])
        row.update(
            external_valid=bool(final), external_probability=probability,
            source_external_probability=source_probability,
            external_change=probability - source_probability,
            external_fold_sd=predictions[final]['fold_sd'][index] if final else None,
            external_fold_probabilities=[v[index] for v in predictions[final]['folds']] if final else [],
            external_score_provenance=predictions[final]['fold_provenance'] if final else [],
            source_external_score_provenance=predictions[source]['fold_provenance'] if source else [],
            operational_route_solved=bool(route_operational),
            operational_solved_external_probability=float(route_operational) * probability,
            source_morgan_similarity=similarities[(source, final)],
            public_overlap_available=public is not None,
            public_label_overlap=bool(matches) if public is not None else None,
            public_task_labeled_overlap=(any(m[TASKS[index]] not in ('', None) for m in matches)
                                        if public is not None else None),
            public_matching_rows=matches,
            output_structure_status=('fixed_final' if row['scientific_endpoint_available']
                                     else 'operational_fallback_scientific_policy_unknown'),
        )
        final_profile, source_profile = profiles.get(final), profiles.get(source)
        row['descriptors'] = final_profile
        row['source_descriptors'] = source_profile
        for name in ('qed', 'molecular_weight'):
            row[name + '_delta'] = (final_profile[name] - source_profile[name]
                                   if final_profile and source_profile else None)
        row['operational_external_joint'] = {
            str(t): bool(route_operational and probability >= t) for t in EXTERNAL_THRESHOLDS}
        row['external_joint'] = {str(t): (row['operational_external_joint'][str(t)]
                                if row['scientific_endpoint_available'] else None)
                                for t in EXTERNAL_THRESHOLDS}
        row['external_joint_bounds'] = {
            str(t): (list(map(int, [row['external_joint'][str(t)]] * 2))
                     if row['scientific_endpoint_available'] else [0, limits(row, 'solved')[1]])
            for t in EXTERNAL_THRESHOLDS}
        row['external_probability_bounds'] = ([probability, probability]
                    if row['scientific_endpoint_available'] else [0., 1.])
        row['solved_external_probability_bounds'] = (
            [row['operational_solved_external_probability']] * 2
            if row['scientific_endpoint_available'] else [0., float(limits(row, 'solved')[1])])
        enriched.append(row)
    return enriched


def distribution(values):
    values = sorted(float(v) for v in values if v is not None)
    if not values:
        return {'n': 0}
    return {'n': len(values), 'mean': statistics.mean(values), 'min': values[0],
            'q25': percentile(values, .25), 'median': percentile(values, .5),
            'q75': percentile(values, .75), 'max': values[-1]}


def filter_pass(row, name):
    profile = row['descriptors']
    if profile is None:
        return False
    if name == 'all':
        return True
    if name == 'no_pains':
        return not profile['pains']
    if name == 'no_pains_brenk':
        return not profile['pains'] and not profile['brenk']
    if name == 'mw_le500':
        return profile['molecular_weight'] <= 500
    if name == 'source_qed_tolerance':
        return profile['qed'] >= row['source_descriptors']['qed'] - .05
    raise ValueError('Unknown structural sensitivity')


def endpoint_bounds(row, threshold, profile_filter='all'):
    lo, hi = row['external_joint_bounds'][str(threshold)]
    if row['scientific_endpoint_available']:
        passed = int(filter_pass(row, profile_filter))
        return lo * passed, hi * passed
    # Unexecuted replacement can have a different external score and descriptors.
    return 0, hi


def rf_quality_bounds(row, profile_filter):
    lo, hi = limits(row, 'joint')
    if profile_filter == 'all':
        return lo, hi
    if row['scientific_endpoint_available']:
        passed = int(filter_pass(row, profile_filter))
        return lo * passed, hi * passed
    return 0, hi


def summarize(rows):
    n = len(rows)
    unknown = [r for r in rows if not r['scientific_endpoint_available']]
    valid = [r for r in rows if r['external_valid']]
    changed = [r for r in valid if r['replaced']]
    joint = {}
    sensitivity = {}
    for threshold in EXTERNAL_THRESHOLDS:
        bounds = [endpoint_bounds(r, threshold) for r in rows]
        joint[str(threshold)] = {
            'operational_count': sum(r['operational_external_joint'][str(threshold)] for r in rows),
            'known_successes': sum(r['external_joint'][str(threshold)] is True for r in rows),
            'scientific_count_bounds': [sum(b[i] for b in bounds) for i in range(2)],
            'scientific_rate_bounds': [sum(b[i] for b in bounds) / n for i in range(2)],
            'unique_known_task_structures': len({(r['task'], r['final_smiles']) for r in rows
                                                if r['external_joint'][str(threshold)] is True}),
        }
        sensitivity[str(threshold)] = {
            name: {
                'operational_count': sum(r['operational_external_joint'][str(threshold)] and
                                         filter_pass(r, name) for r in rows),
                'scientific_count_bounds': [sum(endpoint_bounds(r, threshold, name)[i]
                                              for r in rows) for i in range(2)],
            } for name in FILTERS}
    profiles = [r['descriptors'] for r in valid]
    scaffolds = Counter(p['murcko_scaffold'] for p in profiles if p['murcko_scaffold'])
    return {
        'raw_n': n, 'operational_valid': len(valid), 'invalid': n - len(valid),
        'scientific_endpoint_unknown_slots': len(unknown),
        'unknown_sample_ids': [r['sample_id'] for r in unknown],
        'status_counts': dict(Counter(r['status'] for r in rows)),
        'original_rf_joint_count_bounds': [sum(limits(r, 'joint')[i] for r in rows) for i in range(2)],
        'route_count_bounds': [sum(limits(r, 'solved')[i] for r in rows) for i in range(2)],
        'external_mean_raw_operational': statistics.mean(r['external_probability'] for r in rows),
        'external_mean_raw_scientific_bounds': [
            statistics.mean(r['external_probability_bounds'][i] for r in rows) for i in range(2)],
        'solved_external_mean_raw_operational': statistics.mean(
            r['operational_solved_external_probability'] for r in rows),
        'solved_external_mean_raw_scientific_bounds': [
            statistics.mean(r['solved_external_probability_bounds'][i] for r in rows)
            for i in range(2)],
        'source_external_mean_raw': statistics.mean(r['source_external_probability'] for r in rows),
        'external_source_change_operational': distribution(r['external_change'] for r in rows),
        'probability_increased_slots': sum(r['external_change'] > 1e-8 for r in rows),
        'probability_decreased_slots': sum(r['external_change'] < -1e-8 for r in rows),
        'joint_by_external_threshold': joint, 'quality_joint_sensitivity': sensitivity,
        'original_rf_joint_quality_sensitivity': {
            name: {'scientific_count_bounds': [sum(rf_quality_bounds(r, name)[i]
                                                   for r in rows) for i in range(2)]}
            for name in FILTERS},
        'unique_operational_task_structures': len({(r['task'], r['final_smiles']) for r in valid}),
        'unique_operational_structures': len({r['final_smiles'] for r in valid}),
        'operational_duplicate_task_structure_occupancies': len(valid) -
                            len({(r['task'], r['final_smiles']) for r in valid}),
        'replaced_slots': len(changed),
        'operational_structure_profile': {
            'n_valid_slots': len(valid), 'nonempty_unique_murcko_scaffolds': len(scaffolds),
            'empty_murcko_scaffold_slots': sum(not p['murcko_scaffold'] for p in profiles),
            'largest_same_scaffold_slots': max(scaffolds.values(), default=0),
            'pains_alert_slots': sum(bool(p['pains']) for p in profiles),
            'brenk_alert_slots': sum(bool(p['brenk']) for p in profiles),
            'mw_gt500_slots': sum(p['molecular_weight'] > 500 for p in profiles),
            'qed': distribution(p['qed'] for p in profiles),
            'molecular_weight': distribution(p['molecular_weight'] for p in profiles),
            'heavy_atoms': distribution(p['heavy_atoms'] for p in profiles),
            'source_morgan_similarity': distribution(r['source_morgan_similarity'] for r in valid),
            'changed_source_morgan_similarity': distribution(r['source_morgan_similarity'] for r in changed),
            'qed_delta': distribution(r['qed_delta'] for r in changed),
            'mw_delta': distribution(r['molecular_weight_delta'] for r in changed),
            'changed_same_murcko_scaffold': sum(r['descriptors']['murcko_scaffold'] ==
                         r['source_descriptors']['murcko_scaffold'] for r in changed),
            'interpretation': 'Profiles of fixed operational outputs; unknown slots are '
                              'flagged separately and do not identify the unobserved replacement.',
        },
        'public_overlap_available': all(r['public_overlap_available'] for r in rows),
        'public_overlap_slots': (sum(bool(r['public_label_overlap']) for r in rows)
                                if all(r['public_overlap_available'] for r in rows) else None),
        'public_task_labeled_overlap_slots': (sum(bool(r['public_task_labeled_overlap']) for r in rows)
                                if all(r['public_overlap_available'] for r in rows) else None),
    }


def paired(first, second, runs, threshold, profile_filter='all'):
    if [r['sample_id'] for r in first] != [r['sample_id'] for r in second]:
        raise ValueError('Fixed-output pairing differs')
    effects = {task: [] for task in TASK_INDEX}
    per_run = []
    for spec in runs:
        selected = [(a, b) for a, b in zip(first, second)
                    if a['task'] == spec['task'] and a['replicate'] == spec['replicate']]
        if len(selected) != 64:
            raise ValueError('All64 original slots are needed for each planned training run')
        low = high = wins = losses = unknown = 0
        for a, b in selected:
            al, au = endpoint_bounds(a, threshold, profile_filter)
            bl, bu = endpoint_bounds(b, threshold, profile_filter)
            low += al - bu
            high += au - bl
            if a['scientific_endpoint_available'] and b['scientific_endpoint_available']:
                wins += al == 1 and bl == 0
                losses += bl == 1 and al == 0
            else:
                unknown += 1
        if unknown == 0:
            effects[spec['task']].append(low / 64)
        per_run.append({
            'cohort_id': f"{spec['task']}-rep{spec['replicate']}", 'task': spec['task'],
            'replicate': spec['replicate'], 'raw_n': 64, 'known_wins': wins,
            'known_losses': losses, 'unknown_pairs': unknown,
            'net_count_bounds': [low, high], 'rate_difference_bounds': [low / 64, high / 64],
            'operational_net_count': sum(int(a['operational_external_joint'][str(threshold)] and
                        filter_pass(a, profile_filter)) - int(b['operational_external_joint'][str(threshold)]
                        and filter_pass(b, profile_filter)) for a, b in selected),
        })
        for field, bound in [('external_probability', 'external_probability_bounds'),
                             ('operational_solved_external_probability',
                              'solved_external_probability_bounds')]:
            per_run[-1][field + '_operational_mean_difference'] = statistics.mean(
                a[field] - b[field] for a, b in selected)
            per_run[-1][field + '_scientific_mean_difference_bounds'] = [
                statistics.mean(a[bound][i] - b[bound][1 - i] for a, b in selected)
                for i in range(2)]
    bounds = [sum(r['net_count_bounds'][i] for r in per_run) for i in range(2)]
    result = {'raw_n': 512, 'threshold': threshold, 'profile_filter': profile_filter,
              'runs': per_run, 'net_count_bounds': bounds,
              'rate_difference_bounds': [v / 512 for v in bounds],
              'known_wins': sum(r['known_wins'] for r in per_run),
              'known_losses': sum(r['known_losses'] for r in per_run),
              'unknown_pairs': sum(r['unknown_pairs'] for r in per_run),
              'operational_net_count': sum(r['operational_net_count'] for r in per_run)}
    if result['unknown_pairs'] == 0:
        result['run_level_uncertainty'] = exact_stratified_bootstrap(effects)
        result['positive_runs'] = sum(r['net_count_bounds'][0] > 0 for r in per_run)
        result['negative_runs'] = sum(r['net_count_bounds'][0] < 0 for r in per_run)
    else:
        result['run_level_uncertainty'] = None
        result['uncertainty_note'] = 'Execution bounds retained; no complete-endpoint bootstrap'
    for field, bound in [('external_probability', 'external_probability_bounds'),
                         ('operational_solved_external_probability', 'solved_external_probability_bounds')]:
        result[field + '_operational_mean_difference'] = statistics.mean(a[field] - b[field]
                                                            for a, b in zip(first, second))
        result[field + '_scientific_mean_difference_bounds'] = [
            statistics.mean(a[bound][i] - b[bound][1 - i] for a, b in zip(first, second))
            for i in range(2)]
    result['by_task'] = {
        task: {'raw_n': 256,
               'net_count_bounds': [sum(r['net_count_bounds'][i] for r in per_run
                                       if r['task'] == task) for i in range(2)],
               'unknown_pairs': sum(r['unknown_pairs'] for r in per_run if r['task'] == task)}
        for task in TASK_INDEX}
    return result


def analyze(design, sources, rows, predictions, profiles, similarities, public, audits):
    enriched = enrich(rows, predictions, profiles, similarities, public)
    policies = {p: [r for r in enriched if r['policy'] == p] for p in POLICIES}
    summaries = {p: summarize(selected) for p, selected in policies.items()}
    by_run = {f"{s['task']}-rep{s['replicate']}": {
        p: summarize([r for r in selected if r['task'] == s['task'] and
                      r['replicate'] == s['replicate']])
        for p, selected in policies.items()} for s in design['runs']}
    by_task = {task: {p: summarize([r for r in selected if r['task'] == task])
                     for p, selected in policies.items()} for task in TASK_INDEX}
    contrasts = {}
    for name, left, right in CONTRASTS:
        contrasts[name] = {'first': left, 'second': right,
            'thresholds': {str(t): paired(policies[left], policies[right], design['runs'], t)
                           for t in EXTERNAL_THRESHOLDS}}
    quality_primary = {name: paired(policies['Mix_joint6plus6'], policies['B_joint6plus6'],
                                   design['runs'], .5, name) for name in FILTERS if name != 'all'}
    return {
        'scope': 'raw512_all7_fixed_outputs_cross_architecture_and_structure_sensitivity',
        'raw_denominator': 512, 'source_training_runs': 8, 'policy_slot_rows': len(enriched),
        'model_reference': MODEL_REFERENCE, 'primary_descriptive_threshold': .5,
        'exploratory_thresholds': [.3, .7], 'policies': summaries, 'by_task': by_task,
        'by_training_run': by_run, 'contrasts': contrasts,
        'primary_external_contrast': contrasts[CONTRASTS[0][0]]['thresholds']['0.5'],
        'primary_quality_sensitivities_at_05': quality_primary,
        'public_overlap': ({k: public.get(k) for k in ('source', 'public_rows', 'public_unique_valid',
                                                      'public_invalid_rows')} if public else None),
        'source_unique_task_structures': len({(s['task'], s['canonical_smiles']) for s in sources
                                             if s['valid']}),
        'execution_and_reuse': audits, 'interpretation_limits': LIMITATIONS,
        'rows': enriched, 'fixed_output_note': 'Original rows are copied and enriched only. '
                 'Original RF values, chosen output, Graph endpoints and status remain unchanged.',
    }


def evaluation_design(study, reports, assets, history_paths):
    return {
        'scope': 'raw512_all7_fixed_outputs_external_and_structural_sensitivity',
        'study': str(study), 'source_reports': [str(p) for p in reports],
        'models': MODEL_REFERENCE, 'assets': str(assets),
        'history_paths': [str(p) for p in history_paths], 'primary_descriptive_threshold': .5,
        'exploratory_thresholds': [.3, .7], 'raw_denominator': 512,
        'policy_count': 7, 'training_runs': 8, 'cpu_threads': 2, 'batch_size': 32,
        'missing_only': 'Reuse canonical structure vectors and task-column order per fold; '
            'infer only missing checkpoint/structure pairs and save after each batch.',
        'fixed_outputs': 'All512 source positions x7 policy final outputs exactly; no reselection '
            'or endpoint changes. Source predictions additionally support paired source change.',
        'unknown': 'Operational fallback is recorded. Scientific external score/joint/quality '
            'bounds permit0..1 for unobserved replacement; fallback RF nondegradation does not '
            'imply external-score nondegradation. Invalid sources remain known zeros.',
        'structural_metrics': ['QED', 'MW', 'PAINS', 'BRENK', 'Murcko', 'unique task+canonical',
                               'Morgan radius2/2048 source Tanimoto'],
        'comparisons': CONTRASTS, 'quality_filters': FILTERS,
        'interpretation_limits': LIMITATIONS,
    }


def execute(mode, study, reports, dev, assets, out, chemistry=None, history_paths=None,
            predictor_factory=PublishedPredictor):
    """Prepare/run/analyze; injected test chemistry/predictors never load models."""
    start = time.monotonic()
    chemistry = chemistry or Chemistry()
    design, sources, rows = fixed_inputs(study, reports, chemistry)
    required, requested_pairs = required_structures(sources, rows)
    history_paths = default_history(dev) if history_paths is None else history_paths
    history, provenance, historical_audit = load_history(history_paths, chemistry.canonical)
    out.mkdir(parents=True, exist_ok=True)
    plan = evaluation_design(study, reports, assets, history_paths)
    save(out / 'evaluation_design.json', plan)
    coverage = {}
    for fold in range(10):
        _, current = new_fold(out, fold)
        coverage[str(fold)] = {
            'reused_historical_structures': sum(s in history[fold] for s in required),
            'reused_current_structures': sum(s not in history[fold] and s in current['predictions']
                                            for s in required),
            'missing_structures': [s for s in required if s not in history[fold]
                                   and s not in current['predictions']],
        }
    save(out / 'prediction_plan.json', {
        'tasks': list(TASKS), 'model_reference': MODEL_REFERENCE, 'required_smiles': required,
        'requested_task_canonical_pairs': requested_pairs, 'fold_coverage': coverage,
        'historical_reuse': historical_audit, 'source_reports': [str(p) for p in reports],
    })
    if mode == 'prepare':
        result = {'raw_denominator': 512, 'policy_slot_rows': len(rows),
                  'required_unique_structures': len(required),
                  'missing_fold_structure_predictions': sum(len(v['missing_structures'])
                                                            for v in coverage.values()),
                  'new_model_predictions': 0, 'historical_reuse': historical_audit}
        save(out / 'preparation.json', result)
        return result
    predictions, inference_audit = obtain_predictions(
        required, history, provenance, assets, out, mode == 'run', predictor_factory)
    save(out / 'predictions.json', {'tasks': list(TASKS), 'models': 10,
         'model_reference': MODEL_REFERENCE, **inference_audit, 'predictions': predictions})
    profiles, similarities, descriptor_audit = descriptor_records(required, rows, chemistry, out)
    public, index_path = public_index(dev, assets, out, chemistry)
    result = analyze(design, sources, rows, predictions, profiles, similarities, public, {
        'historical_cache': historical_audit, 'property_inference': inference_audit,
        'descriptor_work': descriptor_audit, 'public_index': index_path,
        'seconds_this_invocation': time.monotonic() - start,
    })
    save(out / 'analysis.json', result)
    with (out / 'slot_outcomes.csv').open('w', encoding='utf-8', newline='') as stream:
        fields = ['policy', 'sample_id', 'task', 'replicate', 'position', 'source_smiles',
                  'final_smiles', 'scientific_endpoint_available', 'status', 'solved', 'joint',
                  'replaced', 'source_score', 'final_score', 'source_qed', 'final_qed',
                  'external_probability', 'source_external_probability', 'external_change',
                  'external_fold_sd', 'source_morgan_similarity', 'qed_delta', 'molecular_weight_delta',
                  'qed', 'molecular_weight', 'heavy_atoms', 'murcko_scaffold', 'pains', 'brenk',
                  'output_structure_status', 'public_label_overlap', 'public_task_labeled_overlap']
        for threshold in EXTERNAL_THRESHOLDS:
            fields.extend(f'external_joint_{threshold}_{suffix}'
                          for suffix in ('known', 'operational', 'lower', 'upper'))
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for row in result['rows']:
            record = dict(row)
            profile = row['descriptors'] or {}
            record.update({k: profile.get(k) for k in
                           ('qed', 'molecular_weight', 'heavy_atoms', 'murcko_scaffold')})
            record.update({k: '|'.join(profile[k]) if k in profile else None
                           for k in ('pains', 'brenk')})
            for threshold in EXTERNAL_THRESHOLDS:
                lower, upper = row['external_joint_bounds'][str(threshold)]
                record.update({f'external_joint_{threshold}_known': row['external_joint'][str(threshold)],
                               f'external_joint_{threshold}_operational':
                                   row['operational_external_joint'][str(threshold)],
                               f'external_joint_{threshold}_lower': lower,
                               f'external_joint_{threshold}_upper': upper})
            writer.writerow(record)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('prepare', 'run', 'analyze', 'self-test'), required=True)
    parser.add_argument('--study', type=Path, default=DEV / 'joint_validation512_v1')
    parser.add_argument('--dev', type=Path, default=DEV)
    parser.add_argument('--assets', type=Path)
    parser.add_argument('--results', type=Path, nargs='+')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if args.mode == 'self-test':
        self_test()
        return
    if not args.results:
        parser.error('--results must provide complete fixed Graph policy reports')
    configure_cpu()
    result = execute(args.mode, args.study, args.results, args.dev,
                     args.assets or args.dev / 'external_property_raw96_v1',
                     args.out or args.study / 'external_fixed_outputs_v1')
    print(json.dumps({key: result[key] for key in (
        ('raw_denominator', 'policy_slot_rows', 'required_unique_structures',
         'missing_fold_structure_predictions', 'new_model_predictions')
        if args.mode == 'prepare' else
        ('raw_denominator', 'policy_slot_rows', 'primary_external_contrast', 'execution_and_reuse'))
    }, indent=2, allow_nan=False))


def self_test():
    """Temporary files, mock chemistry and mock vectors; no real models/imports."""
    import copy
    import tempfile
    from joint_validation512_graph import evaluate_policy

    class MockChemistry:
        version = 'mock-no-RDKit'

        @staticmethod
        def canonical(smiles):
            return {'OCC': 'CCO'}.get(smiles, smiles) if smiles else None

        @staticmethod
        def profile(smiles):
            return {'qed': .60 if smiles == 'CCO' else .57,
                    'molecular_weight': 46. if smiles == 'CCO' else 45.,
                    'heavy_atoms': 3, 'pains': [],
                    'brenk': ['mock_alert'] if smiles == 'CCN' else [],
                    'murcko_scaffold': '', 'components': 1}

        @staticmethod
        def similarity(source, final):
            return None if not source or not final else (1. if source == final else .5)

    def forbidden(*args):
        raise AssertionError('Mock-only test attempted a real or unnecessary model load/search')

    def expect_error(error, function, text):
        try:
            function()
        except error as exc:
            assert text in str(exc), (text, str(exc))
        else:
            raise AssertionError('Expected ' + error.__name__ + ': ' + text)

    checks = []
    with tempfile.TemporaryDirectory(prefix='external512-mock-') as directory:
        base = Path(directory)
        dev, study, assets, out = (base / p for p in ('development', 'study', 'assets', 'out'))
        runs = [{'task': task, 'task_index': index, 'replicate': replicate}
                for replicate in range(4) for index, task in enumerate(TASK_INDEX)]
        design = {'raw_denominator': 512, 'runs': runs,
                  'policies': [{'name': p} for p in POLICIES]}
        sources, rows = [], []
        for spec in runs:
            for position in range(64):
                valid = position != 0
                source = {
                    **spec, 'position': position,
                    'sample_id': f"{spec['task']}-rep{spec['replicate']}-p{position}",
                    'cohort_id': f"{spec['task']}-rep{spec['replicate']}",
                    'canonical_smiles': 'CCO' if valid else None, 'valid': valid,
                    'task_score': .02 if valid else 0., 'qed': .60 if valid else None,
                    'search_seed': 10000 + position, 'graph_verified': valid,
                    'graph_result_path': 'mock/source.json',
                    'graph_result': {'status': 'solved' if valid else 'invalid'},
                }
                sources.append(source)
                for policy in POLICIES:
                    row = evaluate_policy(source, policy, {'scored_candidates': []}, forbidden)
                    if valid and position % 2 and policy == 'Mix_joint6plus6':
                        row.update(final_smiles='CCN', final_score=.8, final_qed=.57,
                                   solved=True, joint=True, replaced=True, status='candidate_solved')
                    if (spec['task'] == 'gsk3b' and spec['replicate'] == 0 and position == 2
                            and policy == 'SynTwins_joint6plus6'):
                        row.update(solved=None, joint=None, scientific_endpoint_available=False,
                                   status='candidate_execution_limitation', endpoint_lower_solved=True,
                                   endpoint_upper_solved=True, endpoint_lower_joint=False,
                                   endpoint_upper_joint=True)
                    rows.append(row)
        original = copy.deepcopy(rows)
        save(study / 'design.json', design)
        save(study / 'sources.json', {'rows': sources})
        report = base / 'graph_report.json'
        save(report, {'rows': rows})
        chemistry = MockChemistry()
        old = dev / 'external_property_raw96_v1/predictions.json'
        folds = [[.8, .2] for _ in range(10)]
        save(old, {'tasks': list(TASKS), 'models': 10,
                   'predictions': {'CCO': {'mean': [.8, .2], 'folds': folds},
                                   'OCC': {'mean': [.8, .2], 'folds': folds}}})
        for fold in range(9):
            save(dev / f'budget_replay256_v1/external_new_fold_{fold:02d}.json',
                 {'fold': fold, 'predictions': {'CCN': [.6, .9]}})
        save(dev / 'external_property_independent256_v1/public_canonical_index.json',
             {'source': 'mock/public.csv', 'public_rows': 1, 'public_unique_valid': 1,
              'public_invalid_rows': 0, 'lookup': {'CCO': [
                  {'public_row': 0, 'gsk3b_activity': '1', 'jnk3_activity': ''}]}})
        history_paths = default_history(dev)
        history, provenance, cache_audit = load_history(history_paths, chemistry.canonical)
        assert cache_audit['historical_unique_structures'] == 2
        assert cache_audit['historical_complete10_structures'] == 1
        assert cache_audit['equivalent_entries'] == 10
        checks.append('canonical and documented legacy tenfold reuse')
        prepared = execute('prepare', study, [report], dev, assets, out,
                           chemistry=chemistry, predictor_factory=forbidden)
        assert prepared['policy_slot_rows'] == 3584 and prepared['missing_fold_structure_predictions'] == 1
        assert read(out / 'prediction_plan.json')['fold_coverage']['9']['missing_structures'] == ['CCN']
        assert not (out / 'analysis.json').exists()
        checks.append('prepare records all512x7 without loading predictors')
        expect_error(RuntimeError, lambda: execute('analyze', study, [report], dev, assets, out,
                     chemistry=chemistry, predictor_factory=forbidden), 'missing1')
        calls = []

        def mock_factory(checkpoint):
            fold = int(checkpoint.parents[1].name.split('_')[1])
            assert fold == 9
            def predict(batch):
                calls.extend((fold, smiles) for smiles in batch)
                assert batch == ['CCN']
                return [[.6, .9]]
            return predict

        result = execute('run', study, [report], dev, assets, out,
                         chemistry=chemistry, predictor_factory=mock_factory)
        assert calls == [(9, 'CCN')]
        assert rows == original and read(report)['rows'] == original
        for before, after in zip(original, result['rows']):
            assert all(after[k] == value for k, value in before.items())
        assert result['raw_denominator'] == 512 and result['policy_slot_rows'] == 3584
        assert len(result['by_training_run']) == 8
        assert all(v['raw_n'] == 512 and v['invalid'] == 8 for v in result['policies'].values())
        known = result['primary_external_contrast']
        assert known['net_count_bounds'] == [128, 128] and known['unknown_pairs'] == 0
        assert known['run_level_uncertainty']['resample_combinations'] == 65536
        assert known['run_level_uncertainty']['interval_95'] == [.25, .25]
        assert len(known['runs']) == 8
        assert abs(known['runs'][0]['external_probability_operational_mean_difference'] + .1) < 1e-8
        with (out / 'slot_outcomes.csv').open(encoding='utf-8', newline='') as stream:
            csv_rows = list(csv.DictReader(stream))
        assert len(csv_rows) == 3584
        csv_unknown = next(r for r in csv_rows if r['scientific_endpoint_available'] == 'False')
        assert csv_unknown['external_joint_0.5_known'] == ''
        assert csv_unknown['external_joint_0.5_operational'] == 'True'
        assert (csv_unknown['external_joint_0.5_lower'], csv_unknown['external_joint_0.5_upper']) == ('0', '1')
        assert any(r['final_smiles'] == 'CCN' and r['brenk'] == 'mock_alert' for r in csv_rows)
        checks.append('fixed outputs preserved; invalid denominator and paired eight-run endpoints')
        unknown = next(r for r in result['rows'] if not r['scientific_endpoint_available'])
        assert unknown['solved'] is None and unknown['joint'] is None
        assert unknown['operational_external_joint']['0.5'] is True
        assert unknown['external_joint']['0.5'] is None
        assert unknown['external_joint_bounds']['0.5'] == [0, 1]
        assert unknown['external_probability_bounds'] == [0., 1.]
        assert unknown['endpoint_lower_solved'] is True
        assert endpoint_bounds(unknown, .5, 'no_pains_brenk') == (0, 1)
        native_contrast = result['contrasts']['B_minus_native']['thresholds']['0.5']
        assert native_contrast['unknown_pairs'] == 1
        assert native_contrast['run_level_uncertainty'] is None
        assert native_contrast['net_count_bounds'] == [0, 1]
        checks.append('execution unknown retains RF/route bounds; external fallback is no lower bound')
        mixture = result['policies']['Mix_joint6plus6']
        assert mixture['operational_structure_profile']['brenk_alert_slots'] == 256
        assert mixture['operational_structure_profile']['empty_murcko_scaffold_slots'] == 504
        assert mixture['unique_operational_task_structures'] == 4
        assert mixture['quality_joint_sensitivity']['0.5']['no_pains_brenk']['scientific_count_bounds'] == [124, 124]
        assert mixture['original_rf_joint_quality_sensitivity']['no_pains_brenk']['scientific_count_bounds'] == [0, 0]
        checks.append('QED/MW/alerts/Murcko/unique/Tanimoto and fullslot sensitivity')
        rerun = execute('analyze', study, [report], dev, assets, out,
                        chemistry=chemistry, predictor_factory=forbidden)
        audit = rerun['execution_and_reuse']['property_inference']
        assert audit['model_structure_predictions_this_invocation'] == 0
        assert audit['loaded_models_this_invocation'] == 0
        assert audit['reused_historical_fold_structure_predictions'] == 19
        assert audit['reused_current_fold_structure_predictions'] == 1
        assert rerun['primary_external_contrast'] == result['primary_external_contrast']
        checks.append('complete-cache resume is inference-free and leaves conclusions identical')
        altered = copy.deepcopy(original)
        altered[1]['source_smiles'] = 'CCN'
        save(report, {'rows': altered})
        expect_error(ValueError, lambda: fixed_inputs(study, [report], chemistry), 'source_smiles')
        altered = copy.deepcopy(original)
        altered[1]['status'] = 'pending_candidate_stage'
        save(report, {'rows': altered})
        expect_error(ValueError, lambda: fixed_inputs(study, [report], chemistry), 'still pending')
        save(report, {'rows': original[:-1]})
        expect_error(ValueError, lambda: fixed_inputs(study, [report], chemistry), 'exactly once')
        save(report, {'rows': original + [original[0]]})
        expect_error(ValueError, lambda: fixed_inputs(study, [report], chemistry), 'exactly once')
        bad = base / 'wrong_tasks.json'
        save(bad, {'tasks': list(reversed(TASKS)), 'models': 10, 'predictions': {}})
        expect_error(ValueError, lambda: load_history([bad], chemistry.canonical), 'task-column')
        conflict = dev / 'budget_replay256_v1/conflict.json'
        save(conflict, {'fold': 0, 'predictions': {'CCO': [.7, .2]}})
        expect_error(ValueError, lambda: load_history([old, conflict], chemistry.canonical), 'Conflicting')
        checks.append('ordinary source alignment, pending/coverage and cache-conflict rejection')

        required = ['mock' + str(i).zfill(3) for i in range(65)]
        restart_out = base / 'restart'
        empty = [{} for _ in range(10)]
        interrupted_calls = []

        def interrupting_factory(checkpoint):
            assert checkpoint.parents[1].name == 'fold_0'
            batches = 0
            def predict(batch):
                nonlocal batches
                batches += 1
                if batches == 2:
                    raise KeyboardInterrupt('mock batch interrupt')
                interrupted_calls.extend((0, s) for s in batch)
                return [[.4, .6] for _ in batch]
            return predict

        expect_error(KeyboardInterrupt, lambda: obtain_predictions(required, empty, empty, assets,
                     restart_out, True, interrupting_factory), 'mock batch interrupt')
        assert len(read(restart_out / 'new_fold_00_predictions.json')['predictions']) == 32
        resumed_calls = []

        def resumed_factory(checkpoint):
            fold = int(checkpoint.parents[1].name.split('_')[1])
            def predict(batch):
                resumed_calls.extend((fold, s) for s in batch)
                return [[.4, .6] for _ in batch]
            return predict

        resumed_predictions, resume_audit = obtain_predictions(
            required, empty, empty, assets, restart_out, True, resumed_factory)
        assert len(resumed_predictions) == 65
        assert not set(interrupted_calls) & set(resumed_calls)
        assert len(interrupted_calls + resumed_calls) == 650
        assert len(set(interrupted_calls + resumed_calls)) == 650
        assert resume_audit['model_structure_predictions_this_invocation'] == 618
        assert resume_audit['reused_current_fold_structure_predictions'] == 32
        _, no_work = obtain_predictions(required, empty, empty, assets, restart_out, False, forbidden)
        assert no_work['loaded_models_this_invocation'] == 0
        assert no_work['model_structure_predictions_this_invocation'] == 0
        assert no_work['reused_current_fold_structure_predictions'] == 650
        checks.append('mid-fold interrupt saves32; restart completes650 pairs without repeats')
    print(json.dumps({'self_test': 'passed', 'checks': checks, 'real_model_calls': 0,
                      'real_search_calls': 0, 'mock_slot_rows': 3584}, indent=2))


if __name__ == '__main__':
    main()

