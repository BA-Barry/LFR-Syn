"""Core extracted from retained v2_unified_raw8_candidates.py; see code_origins.json."""
from collections import Counter
import time
from rdkit import Chem
from . import common
NATIVE = {}


def find_roots(source, inventory, forward, retros):
    started = time.monotonic()
    molecule = Chem.MolFromSmiles(source['canonical_smiles'])
    seen, roots, costs = set(), [], Counter()
    for ref, retro in retros:
        if forward[ref].GetNumReactantTemplates() != 2:
            continue
        costs['retro_template_calls'] += 1
        try:
            pairs = retro.RunReactants((molecule,))
        except Exception:
            costs['retro_errors'] += 1
            continue
        for pair in pairs:
            precursors = [common.canon(Chem.Mol(m)) for m in pair]
            if len(precursors) != 2 or not all(precursors):
                continue
            key = (ref, tuple(precursors))
            if key in seen:
                continue
            seen.add(key)
            costs['source_replay_forward_calls'] += 1
            try:
                recovered = any(source['canonical_smiles'] in {common.canon(Chem.Mol(p)) for p in ps}
                    for ps in forward[ref].RunReactants(tuple(Chem.MolFromSmiles(s) for s in precursors)))
            except Exception:
                costs['source_replay_errors'] += 1
                recovered = False
            if recovered:
                missing = [i for i,s in enumerate(precursors) if Chem.MolToInchiKey(Chem.MolFromSmiles(s)) not in NATIVE['stock_keys']]
                roots.append({'rxn_ref': ref, 'precursors': precursors, 'missing_sides': missing})
    groups = {1: [], 2: []}
    for count in (1, 2):
        used_refs = set()
        for r in sorted(roots, key=lambda r: (r['rxn_ref'], tuple(r['precursors']))):
            if len(r['missing_sides']) == count and r['rxn_ref'] not in used_refs:
                groups[count].append(r)
                used_refs.add(r['rxn_ref'])
    selected = groups[1][:2] + groups[2][:2]
    extras = groups[1][2:] + groups[2][2:]
    selected.extend(extras[:4-len(selected)])
    selected.sort(key=lambda r: (len(r['missing_sides']), r['rxn_ref'], tuple(r['precursors'])))
    counts = Counter(len(r['missing_sides']) for r in roots)
    return selected, dict(costs, root_seconds=time.monotonic()-started,
        exact_replay_root_counts=dict(counts), grouped_root_counts={str(k):len(v) for k,v in groups.items()},
        selected_root_counts=dict(Counter(len(r['missing_sides']) for r in selected)),
        all_exact_replay_roots=roots)
