"""Plot saved prospective rates and intervals; perform no statistical inference."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def saved_numbers(analysis_path):
    analysis = json.loads(Path(analysis_path).read_text(encoding='utf-8-sig'))
    if analysis['independent_optimizer_runs'] != 16 or analysis['raw_source_positions'] != 1024:
        raise ValueError('Figure2 requires the completed sixteen-run study')
    pairs = {}
    for task in ('gsk3b', 'jnk3'):
        pairs[task] = []
        for replicate in range(8):
            cohort = 'graphga-qualconfirm-%s-rep%d' % (task, replicate)
            arm = analysis['by_cohort'][cohort]
            values = []
            for policy in ('B_joint6plus6', 'B_fixedtask6plus6'):
                endpoint = arm[policy]['external_joint']
                if arm[policy]['raw_n'] != 64 or endpoint['unknown_slots']:
                    raise ValueError('Each plotted paired run must retain64known raw positions')
                values.append({'policy': policy, 'known_successes': endpoint['known_successes'],
                               'raw_denominator': 64, 'rate_percent': 100 * endpoint['raw_rate_bounds'][0]})
            pairs[task].append({'cohort_id': cohort, 'replicate': replicate, 'values': values})
    requested = [
        ('B qualification (primary)', analysis['primary']),
        ('C qualification (exploratory)', next(c for c in analysis['exploratory_secondary_contrasts']
          if c['endpoint'] == 'external_joint' and c['terms'] == [
              {'policy': 'C_fixedtask6plus6', 'weight': 1}, {'policy': 'C_joint6plus6', 'weight': -1}])),
        ('Fixed B − fixed C (exploratory)', next(c for c in analysis['exploratory_secondary_contrasts']
          if c['endpoint'] == 'external_joint' and c['terms'] == [
              {'policy': 'B_fixedtask6plus6', 'weight': 1}, {'policy': 'C_fixedtask6plus6', 'weight': -1}])),
    ]
    contrasts = []
    for label, contrast in requested:
        if contrast['technical_unknown_pairs'] or len(contrast['runs']) != 16:
            raise ValueError('Saved contrast is not the complete paired sixteen-run contrast')
        contrasts.append({'label': label, 'terms': contrast['terms'], 'scope': contrast['scope'],
                          'effect_percentage_points': contrast['welch']['task_equal_mean_effect_pp'],
                          'interval_95_percentage_points': contrast['welch']['interval_95_pp'],
                          'interval_method': contrast['welch']['method'], 'optimizer_runs': 16})
    return {'analysis_source': str(Path(analysis_path).resolve()), 'endpoint': 'external_joint',
        'endpoint_definition': 'Unrounded mean intended-task probability across10saved checkpoints >=0.5 and saved audited route; no guiding-RF conjunct.',
        'source_arrays': 'analysis.by_cohort for rates; saved analysis.primary / exploratory_secondary_contrasts Welch entries for point and interval.',
        'paired_runs_by_task': pairs, 'contrasts': contrasts,
        'plotting_only': 'All rates and intervals are read from the completed analysis. Replicate-specific horizontal jitter separates overlapping points, with unchanged y values. No new statistic, confidence interval, model prediction or chemistry is computed.'}


def plot(numbers, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
        'axes.titlesize': 9, 'axes.labelsize': 8, 'xtick.labelsize': 7.5,
        'ytick.labelsize': 7.5, 'svg.fonttype': 'none', 'savefig.facecolor': 'white',
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': .65, 'xtick.major.width': .65, 'ytick.major.width': .65})
    fig = plt.figure(figsize=(7.25, 3.10), layout='constrained')
    grid = fig.add_gridspec(1, 3, width_ratios=(1.0, 1.0, 1.95), wspace=.16)
    blue, grey = '#176B91', '#78828A'
    for index, (task, title) in enumerate((('gsk3b', 'GSK3B'), ('jnk3', 'JNK3'))):
        ax = fig.add_subplot(grid[0, index])
        entries = numbers['paired_runs_by_task'][task]
        for replicate, entry in enumerate(entries):
            # Deterministic display offsets, not new observations or inference.
            offset = (replicate - 3.5) * .017
            rates = [value['rate_percent'] for value in entry['values']]
            xs = [offset, 1 + offset]
            ax.plot(xs, rates, color='#BBC4CA', linewidth=.9, zorder=2)
            ax.scatter(xs[0], rates[0], s=17, facecolor='white', edgecolor=grey, linewidth=.8, zorder=3)
            ax.scatter(xs[1], rates[1], s=17, facecolor=blue, edgecolor=blue, linewidth=.6, zorder=3)
        ax.set(xlim=(-.23, 1.23), ylim=(-4, 106), xticks=(0, 1),
               xticklabels=('Source-\npreserving', 'Fixed task'), yticks=(0, 25, 50, 75, 100))
        ax.set_title(title, loc='center', pad=8, weight='bold')
        ax.grid(axis='y', color='#E6EAED', linewidth=.65, zorder=0)
        ax.tick_params(axis='x', length=0, pad=7)
        ax.text(.5, -.26, '8 paired runs × 64 positions', ha='center', va='top',
                transform=ax.transAxes, fontsize=7, color='#48545D')
        ax.text(-.13, 1.06, chr(ord('a') + index), transform=ax.transAxes,
                fontsize=11, fontweight='bold')
        if index == 0:
            ax.set_ylabel('Predicted qualification + audited route (%)', labelpad=7)
        else:
            ax.set_yticklabels([])
    ax = fig.add_subplot(grid[0, 2])
    for y, item in zip((2, 1, 0), numbers['contrasts']):
        point = item['effect_percentage_points']
        lower, upper = item['interval_95_percentage_points']
        color = blue if y == 2 else grey
        ax.errorbar(point, y, xerr=[[point - lower], [upper - point]], fmt='o',
                    color=color, markersize=4.2, capsize=3, elinewidth=1.4,
                    markerfacecolor=color if y == 2 else 'white', markeredgewidth=1.1, zorder=3)
        ax.text(point, y + .20, '%.1f [%.1f, %.1f]' % (point, lower, upper),
                ha='center', va='bottom', fontsize=7.2, color=color)
    ax.axvline(0, color='#525E66', linewidth=.8, zorder=1)
    ax.set(xlim=(-15, 100), ylim=(-.55, 2.7), yticks=(2, 1, 0),
           yticklabels=('B qualification\n(primary)', 'C qualification\n(exploratory)',
                        'Fixed B − fixed C\n(exploratory)'), xticks=(0, 25, 50, 75, 100))
    ax.tick_params(axis='y', length=0, pad=7)
    ax.grid(axis='x', color='#E6EAED', linewidth=.65, zorder=0)
    ax.set_title('Paired-run effects (95% CI)', pad=8, fontsize=8.7, weight='bold')
    ax.set_xlabel('Difference (percentage points)', labelpad=7)
    ax.text(.5, -.26, '16 optimizer runs; secondary intervals unadjusted',
            ha='center', va='top', transform=ax.transAxes, fontsize=7, color='#48545D')
    ax.text(-.08, 1.06, 'c', transform=ax.transAxes, fontsize=11, fontweight='bold')
    ax.spines['left'].set_visible(False)
    fig.savefig(output / 'figure2_independent_qualification.png', dpi=600, bbox_inches='tight', pad_inches=.05)
    fig.savefig(output / 'figure2_independent_qualification.svg', bbox_inches='tight', pad_inches=.05)
    plt.close(fig)
    (output / 'figure2_source_numbers.json').write_text(json.dumps(numbers, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    numbers = saved_numbers(args.analysis)
    plot(numbers, args.output)
    print(json.dumps({'output': str(args.output.resolve()), 'paired_runs': 16,
                      'run_raw_denominator': 64, 'intervals': numbers['contrasts']}, indent=2, ensure_ascii=True))


if __name__ == '__main__':
    main()
