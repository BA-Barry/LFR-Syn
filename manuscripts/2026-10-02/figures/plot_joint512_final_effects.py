"""Draw complete raw512 policy results and fixed-output external sensitivity.

Every raw slot is retained. Intervals use the saved conditional paired-run
bootstrap; native execution bounds are distinct from statistical intervals.
There is no molecule sampling, model prediction or candidate reselection.
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


POLICIES = ('B_route12', 'B_joint12', 'B_joint6plus6', 'Mix_joint6plus6',
            'C_joint6plus6', 'A_joint6plus6', 'SynTwins_joint6plus6')
LABELS = ('Original', 'B route / 12', 'B joint / 12', 'B joint / 6+6',
          'Mix joint / 6+6', 'C joint / 6+6', 'A joint / 6+6', 'SynTwins / 6+6')
COLORS = ('#BBBBBB', '#82B1CC', '#4D93BB', '#0072B2', '#009E73', '#CCA24A', '#808080', '#8A7EA4')
COMPARISONS = (
    ('primary_Mix_minus_B_split', 'Mix minus B (primary)'),
    ('attribution_B_minus_C', 'B minus C'),
    ('trigger_joint_minus_route', 'Joint minus route trigger'),
    ('allocation_split_minus_twelve', '6+6 minus first-candidate 12'),
    ('B_minus_Morgan', 'B minus A'),
)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--external', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    analysis, external = read(args.analysis), read(args.external)
    if analysis['raw_denominator'] != 512 or analysis['source_training_runs'] != 8:
        raise ValueError('Need the full eight-run 512-slot analysis')
    if external['raw_denominator'] != 512 or external['policy_slot_rows'] != 3584:
        raise ValueError('Need all seven complete fixed-output policies')
    if any('pending_candidate_stage' in analysis['policies'][p]['statuses'] for p in POLICIES):
        raise ValueError('Final figure cannot contain an unmeasured candidate stage')
    originals = [r for r in external['rows'] if r['policy'] == 'B_route12']
    if len(originals) != 512 or len({r['sample_id'] for r in originals}) != 512:
        raise ValueError('Original source coverage is not complete')
    external_original = sum(bool(r['source_solved'] and r['source_external_probability'] >= .5) for r in originals)
    original_joint = analysis['source_counts']['original_joint']
    rf_bounds = [[original_joint, original_joint]] + [analysis['policies'][p]['joint_bounds'] for p in POLICIES]
    external_bounds = [[external_original, external_original]] + [external['policies'][p]['joint_by_external_threshold']['0.5']['scientific_count_bounds'] for p in POLICIES]
    args.output.mkdir(parents=True, exist_ok=True)
    counts = []
    for policy, label, rf, ext in zip(('Original',) + POLICIES, LABELS, rf_bounds, external_bounds):
        counts.append({'policy': policy, 'label': label, 'raw_n': 512,
                       'rf_joint_lower': rf[0], 'rf_joint_upper': rf[1],
                       'external05_joint_lower': ext[0], 'external05_joint_upper': ext[1]})
    write_csv(args.output / 'joint512_final_policy_counts.csv', counts)
    rf_contrasts = {r['name']: r for r in [analysis['primary_contrast']] + analysis['secondary_contrasts']}
    effects = []
    for name, label in COMPARISONS:
        for endpoint, record in [('Optimized RF', rf_contrasts[name]),
                                 ('Chemprop >= 0.5', external['contrasts'][name]['thresholds']['0.5'])]:
            if record['unknown_pairs']:
                raise ValueError('Selected complete-policy contrasts have unknown pairs')
            uncertainty = record['run_level_uncertainty']
            effects.append({'comparison': name, 'label': label, 'endpoint': endpoint, 'raw_n': 512,
                            'difference_percentage_points': 100 * uncertainty['effect'],
                            'interval_lower_percentage_points': 100 * uncertainty['interval_95'][0],
                            'interval_upper_percentage_points': 100 * uncertainty['interval_95'][1],
                            'positive_runs': record['positive_runs'], 'negative_runs': record['negative_runs']})
    write_csv(args.output / 'joint512_final_paired_effects.csv', effects)
    run_rows = []
    for name, label in COMPARISONS:
        for endpoint, record, field in [('Optimized RF', rf_contrasts[name], 'net_joint_bounds'),
            ('Chemprop >= 0.5', external['contrasts'][name]['thresholds']['0.5'], 'net_count_bounds')]:
            for run in record['runs']:
                run_rows.append({'comparison': name, 'endpoint': endpoint, 'cohort_id': run['cohort_id'],
                                 'raw_n': 64, 'net_lower_slots': run[field][0], 'net_upper_slots': run[field][1]})
    write_csv(args.output / 'joint512_final_run_differences.csv', run_rows)

    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
        'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 8,
        'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .6,
        'legend.frameon': False, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(7.086614173, 5.708661417))  # 180 x 145 mm
    top = [fig.add_axes([.155, .635, .29, .295]), fig.add_axes([.655, .635, .29, .295])]
    for ax, bounds, title in zip(top, (rf_bounds, external_bounds),
            ('RF task cutoff + audited route', 'Chemprop ≥ 0.5 + audited route')):
        for y, ((lower, upper), color) in enumerate(zip(bounds, COLORS)):
            ax.barh(y, 100 * lower / 512, height=.55, color=color)
            if upper > lower:
                ax.barh(y, 100 * (upper - lower) / 512, left=100 * lower / 512,
                        height=.55, color=color, alpha=.3, edgecolor=color, linewidth=.6)
            label = str(lower) if lower == upper else f'{lower}–{upper}'
            ax.text(100 * upper / 512 + 1.5, y, label, va='center', fontsize=6.5)
        ax.set(yticks=list(range(8)), xlim=(0, 82), ylim=(7.65, -.65), xticks=[0, 20, 40, 60, 80])
        ax.set_yticklabels(LABELS)
        ax.tick_params(axis='y', length=0, pad=5)
        ax.spines['left'].set_visible(False)
        ax.set_xlabel('Joint outputs / 512 raw slots (%)', labelpad=4)
        ax.set_title(title, pad=10)
    fig.text(.02, .963, 'a', fontsize=9, fontweight='bold')
    fig.text(.515, .963, 'b', fontsize=9, fontweight='bold')
    fig.text(.155, .543, 'Counts shown at bar ends; pale extensions are the eight native execution-unknown slots.', fontsize=6.5)
    fig.text(.02, .507, 'c', fontsize=9, fontweight='bold')
    ax = fig.add_axes([.335, .145, .60, .295])
    ax.axvline(0, color='#A0A0A0', linewidth=.6)
    for index, (name, label) in enumerate(COMPARISONS):
        for endpoint_index, (endpoint, color, marker) in enumerate([
                ('Optimized RF', '#0072B2', 'o'), ('Chemprop >= 0.5', '#D55E00', 's')]):
            r = next(row for row in effects if row['comparison'] == name and row['endpoint'] == endpoint)
            y = index + (-.13 if endpoint_index == 0 else .13)
            lo, hi = r['interval_lower_percentage_points'], r['interval_upper_percentage_points']
            ax.plot([lo, hi], [y, y], color=color, linewidth=1)
            ax.scatter(r['difference_percentage_points'], y, color=color, marker=marker, s=20,
                       edgecolors='none', label=endpoint if index == 0 else None, zorder=3)
    ax.set(yticks=list(range(5)), ylim=(4.6, -.6), xlim=(-9, 28), xticks=[-5, 0, 5, 10, 15, 20, 25])
    ax.set_yticklabels([label for _, label in COMPARISONS])
    ax.tick_params(axis='y', length=0, pad=7)
    ax.spines['left'].set_visible(False)
    ax.set_title('Paired effects across eight complete training runs', pad=10)
    ax.set_xlabel('Difference in joint yield (percentage points)', labelpad=6)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, [label.replace('>=', '≥') for label in labels], loc='center right',
               bbox_to_anchor=(.95, .5), ncol=2, fontsize=6.5, handletextpad=.5)
    fig.text(.155, .055, 'Lines: conditional 95% empirical bootstrap within each fixed task (4 runs per task).', fontsize=6.5)
    fig.text(.155, .033, 'All 512 raw slots retained; predictor sensitivity does not establish experimental activity.', fontsize=6.5)
    base = args.output / 'joint512_final_effects'
    fig.savefig(base.with_suffix('.pdf'))
    fig.savefig(base.with_suffix('.svg'))
    fig.savefig(base.with_suffix('.tiff'), dpi=1200, pil_kwargs={'compression': 'tiff_lzw'})
    fig.savefig(base.with_suffix('.png'), dpi=300)
    plt.close(fig)
    print(json.dumps({'figure': str(base), 'raw_slots': 512, 'complete_policies': 7,
                      'original_external_joint05': external_original}, indent=2))


if __name__ == '__main__':
    main()
