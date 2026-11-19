"""Plot completed source and candidate qualification data, not final route effects.

The input contains all eight independent terminal runs, with 64 raw generation
slots per run. Every mark is an observed count or a deterministic paired
difference. No sampling, smoothing, model inference or outcome selection.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


ARMS = ('morgan', 'permuted_contribution', 'task_contribution', 'mix')
LABELS = ('A: Morgan', 'C: Permuted', 'B: Contribution', 'Mix: B + A')
COLORS = ('#808080', '#CC9F4A', '#0072B2', '#009E73')
COMPARISONS = (
    ('B minus A', 'task_contribution', 'morgan'),
    ('B minus C', 'task_contribution', 'permuted_contribution'),
    ('Mix minus B', 'mix', 'task_contribution'),
)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_csv(path, records):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original', type=Path, required=True)
    parser.add_argument('--candidates', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    original, candidates = read(args.original), read(args.candidates)
    runs = [f'{task}-rep{rep}' for task in ('gsk3b', 'jnk3') for rep in range(4)]
    raw = original['all']['raw']
    if raw != 512 or candidates['raw_denominator'] != raw:
        raise ValueError('Expected the complete 512-slot source cohort')
    if set(original['by_run']) != set(runs) or set(candidates['by_run']) != set(runs):
        raise ValueError('Expected all eight prespecified terminal runs')
    if any(original['by_run'][run]['raw'] != 64 for run in runs):
        raise ValueError('A terminal run has a different raw denominator')
    for arm in ARMS:
        if sum(candidates['by_run'][run][arm]['eligible_slots'] for run in runs) != candidates['arms'][arm]['candidate_eligible_slots']:
            raise ValueError(f'Inconsistent aggregation for {arm}')
    args.output.mkdir(parents=True, exist_ok=True)
    source_records = []
    for run in runs:
        r = original['by_run'][run]
        source_records.append({'run': run, 'task': run.split('-')[0], 'raw_slots': r['raw'],
                               'valid': r['valid'], 'RF_qualified': r['RF_qualified'],
                               'audited_route': r['audited_route'], 'original_joint': r['original_joint'],
                               **r['strata']})
    write_csv(args.output / 'joint512_source_counts.csv', source_records)
    qualification_records = []
    difference_records = []
    for run in runs:
        for arm, label in zip(ARMS, LABELS):
            r = candidates['by_run'][run][arm]
            qualification_records.append({'run': run, 'task': run.split('-')[0], 'arm': arm,
                                           'label': label, 'raw_slots': 64,
                                           'eligible_slots': r['eligible_slots'],
                                           'at_least_two_eligible': r['at_least_two_eligible'],
                                           'eligible_per_raw_percent': 100 * r['eligible_slots'] / 64})
        for label, first, second in COMPARISONS:
            difference = candidates['by_run'][run][first]['eligible_slots'] - candidates['by_run'][run][second]['eligible_slots']
            difference_records.append({'run': run, 'task': run.split('-')[0], 'comparison': label,
                                       'raw_slots': 64, 'difference_slots': difference,
                                       'difference_percentage_points': 100 * difference / 64})
    write_csv(args.output / 'joint512_candidate_qualification.csv', qualification_records)
    write_csv(args.output / 'joint512_candidate_paired_differences.csv', difference_records)

    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
        'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7,
        'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
        'axes.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False,
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'legend.frameon': False,
        'lines.linewidth': 0.8,
    })
    # Physical dimensions are fixed; use final PDF for text/collision QA.
    fig = plt.figure(figsize=(7.086614173, 5.393700787))  # 180 x 137 mm
    a = fig.add_axes([0.095, 0.685, 0.28, 0.225])
    b = fig.add_axes([0.555, 0.685, 0.41, 0.225])
    fig.text(0.02, 0.96, 'a', fontsize=9, fontweight='bold')
    fig.text(0.48, 0.96, 'b', fontsize=9, fontweight='bold')
    fig.text(0.02, 0.57, 'c', fontsize=9, fontweight='bold')
    a.set_title('Original source outcomes', pad=11)
    cells = (
        (0, 0, 'high_task_route_solved', '#D5E9F4'),
        (1, 0, 'high_task_route_unresolved', '#F2F2F2'),
        (0, 1, 'low_task_route_solved', '#F2F2F2'),
        (1, 1, 'low_task_route_unresolved', '#F2F2F2'),
    )
    for x, y, key, color in cells:
        a.add_patch(Rectangle((x, y), 1, 1, facecolor=color, edgecolor='white', linewidth=2))
        count = original['all']['strata'][key]
        a.text(x + 0.5, y + 0.42, str(count), ha='center', va='center', fontsize=11)
        a.text(x + 0.5, y + 0.7, f'{100 * count / raw:.1f}%', ha='center', va='center', fontsize=6.5)
    a.set(xlim=(0, 2), ylim=(2, 0), xticks=[0.5, 1.5], yticks=[0.5, 1.5])
    a.set_xticklabels(['Route solved', 'Unresolved'])
    a.set_yticklabels(['RF qualified', 'Below cutoff'])
    a.tick_params(length=0, pad=4)
    for spine in a.spines.values():
        spine.set_visible(False)
    fig.text(0.095, 0.63, f'{raw} raw slots; {original["all"]["invalid"]} invalid slots outside matrix', fontsize=6.5)
    fig.text(0.095, 0.607, f'{original["all"]["original_joint"]}/{raw} original joint successes (10.2%)', fontsize=6.5)

    b.set_title('Slots with at least one eligible candidate', pad=11)
    for y, arm, label, color in zip(range(4), ARMS, LABELS, COLORS):
        count = candidates['arms'][arm]['candidate_eligible_slots']
        value = 100 * count / raw
        b.barh(y, value, height=0.55, color=color)
        b.text(value + 1.8, y, f'{count}/{raw}', ha='left', va='center', fontsize=7)
    b.set(yticks=list(range(4)), xlim=(0, 76), ylim=(3.65, -0.65), xticks=[0, 20, 40, 60])
    b.set_yticklabels(LABELS)
    b.set_xlabel('Eligible candidate slots / raw slots (%)', labelpad=4)
    b.tick_params(axis='y', length=0, pad=5)
    b.spines['left'].set_visible(False)

    offsets = [0.11, 0.41, 0.71]
    for index, ((title, first, second), left) in enumerate(zip(COMPARISONS, offsets)):
        ax = fig.add_axes([left, 0.13, 0.255, 0.365])
        ax.axvline(0, color='#999999', linewidth=0.65, zorder=1)
        ax.axhspan(3.5, 7.5, color='#F3F3F3', zorder=0)
        for y, run in enumerate(runs):
            difference = candidates['by_run'][run][first]['eligible_slots'] - candidates['by_run'][run][second]['eligible_slots']
            ax.scatter(100 * difference / 64, y, s=20,
                       color='#0072B2' if y < 4 else '#D55E00',
                       marker='o' if y < 4 else 's', edgecolors='none', zorder=2)
        ax.set(xlim=(-18, 36), ylim=(7.55, -0.55), xticks=[-15, 0, 15, 30], yticks=list(range(8)))
        ax.set_title(title, pad=8)
        ax.set_xlabel('Difference (percentage points)', labelpad=5)
        ax.set_yticklabels([f'{"GSK3B" if y < 4 else "JNK3"} r{y % 4 + 1}' for y in range(8)] if index == 0 else [])
        ax.tick_params(axis='y', length=0, pad=5)
        ax.spines['left'].set_visible(False)
    fig.text(0.11, 0.035, 'Each point is a paired difference from one full terminal run (64 raw slots); 4 runs per task.', fontsize=6.5)
    base = args.output / 'joint512_candidate_diagnostics'
    fig.savefig(base.with_suffix('.pdf'))
    fig.savefig(base.with_suffix('.svg'))
    fig.savefig(base.with_suffix('.tiff'), dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    fig.savefig(base.with_suffix('.png'), dpi=300)
    plt.close(fig)
    print(json.dumps({'figure': str(base), 'raw_slots': raw, 'independent_runs': len(runs),
                      'endpoint': 'candidate_qualification_before_candidate_route_audit'}, indent=2))


if __name__ == '__main__':
    main()
