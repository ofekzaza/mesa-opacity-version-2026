import glob
import os
import re
import sys

import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np

import utils

BASE_DIR = '.'

DIRS = [
    "10m-ours_reduction", "10m-ours_reduction_post",
    "10m-supereduction_a=2_reduction", "10m-supereduction_a=2_reduction_post",
    "20m-ours_reduction", "20m-ours_reduction_post",
    "20m-supereduction_a=2_reduction", "20m-supereduction_a=2_reduction_post",
    "30m-ours_reduction", "30m-ours_reduction_post",
    "30m-supereduction_a=2_reduction", "30m-supereduction_a=2_reduction_post",
    "30m-normal_reduction", "30m-normal_reduction_post",
    "30m-mlt++_reduction", "30m-mlt++_reduction_post",
    "40m-ours_reduction", "40m-ours_reduction_post",
    "40m-supereduction_a=2_reduction", "40m-supereduction_a=2_reduction_post",
    "50m-ours_reduction", "50m-ours_reduction_post",
    "50m-supereduction_a=2_reduction", "50m-supereduction_a=2_reduction_post",
    "60m-ours_reduction", "60m-ours_reduction_post",
    "60m-supereduction_a=2_reduction", "60m-supereduction_a=2_reduction_post",
    "70m-ours_reduction", "70m-ours_reduction_post",
    "70m-supereduction_a=2_reduction", "70m-supereduction_a=2_reduction_post",
    "80m-ours_reduction", "80m-ours_reduction_post",
    "80m-supereduction_a=2_reduction", "80m-supereduction_a=2_reduction_post",
    "90m-ours_reduction", "90m-ours_reduction_post",
    "90m-supereduction_a=2_reduction", "90m-supereduction_a=2_reduction_post",
    "100m-ours_reduction", "100m-ours_reduction_post",
    "100m-supereduction_a=2_reduction", "100m-supereduction_a=2_reduction_post",
]


def scan_neg_profiles(log_dir):
    """Return list of (profile_num, model_number, neg_zones) for profiles with negative gradT."""
    prof_files = sorted(
        glob.glob(os.path.join(log_dir, 'profile*.data')),
        key=lambda x: int(x.split('profile')[-1].split('.data')[0])
    )
    results = []
    for fpath in prof_files:
        try:
            prof = mr.MesaData(fpath)
        except Exception:
            continue
        if not hasattr(prof, 'gradT'):
            continue
        gradt = np.array(prof.gradT, dtype=float)
        neg_mask = gradt < 0
        if not neg_mask.any():
            continue
        profile_num = int(os.path.basename(fpath).split('profile')[1].split('.data')[0])
        zones = np.array(prof.zone, dtype=int)[neg_mask]
        model_number = int(prof.header('model_number'))
        results.append((profile_num, model_number, zones, fpath))
    return results


def plot_dir():
    out_dir = 'plots_neg_gradT'
    os.makedirs(out_dir, exist_ok=True)

    for d in DIRS:
        log_dir = os.path.join(BASE_DIR, d, 'LOGS')
        if not os.path.isdir(log_dir):
            continue

        negs = scan_neg_profiles(log_dir)
        if not negs:
            continue

        # Create per-run subdirectory
        run_dir = os.path.join(out_dir, d)
        os.makedirs(run_dir, exist_ok=True)

        print(f"{d}: {len(negs)} profiles with negative gradT", file=sys.stderr)

        for profile_num, model_number, neg_zones, fpath in negs:
            prof = mr.MesaData(fpath)
            logT = np.array(prof.logT, dtype=float)
            logP = np.array(prof.logP, dtype=float)
            gradT = np.array(prof.gradT, dtype=float)
            zones = np.array(prof.zone, dtype=int)

            fig, ax = plt.subplots(figsize=(8, 6))

            # Plot full profile as one continuous line
            ax.plot(logP, logT, color='#1f77b4', linewidth=1.0, label='full profile', zorder=2)

            # Find contiguous negative blocks and overlay them
            neg_mask = gradT < 0
            diff = np.diff(np.concatenate(([False], neg_mask, [False])).astype(int))
            starts = np.where(diff == 1)[0]
            ends = np.where(diff == -1)[0]

            for s, e in zip(starts, ends):
                # Extend one point on each side for visual continuity
                seg_start = max(0, s - 1)
                seg_end = min(len(logP), e + 1)
                ax.plot(logP[seg_start:seg_end], logT[seg_start:seg_end],
                        color='#d62728', linewidth=2.0, zorder=3)

            ax.plot([], [], color='#d62728', linewidth=2.0, label='gradT < 0')

            ax.set_xlabel(utils.pretty_axis_name('logP'))
            ax.set_ylabel(utils.pretty_axis_name('logT'))
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)

            neg_str = ', '.join(str(z) for z in neg_zones[:10])
            if len(neg_zones) > 10:
                neg_str += f', ... ({len(neg_zones)} total)'
            ax.set_title(f'{utils.parse_name_for_plots(d)}  Model {model_number}  (profile {profile_num})\n'
                        f'Neg gradT at layers: {neg_str}', fontsize=9)

            plt.tight_layout()
            plt.savefig(os.path.join(run_dir, f'model{model_number}_prof{profile_num}.png'), dpi=200)
            plt.close()

    # Summary plot: overlay all negative regions
    print("\nGenerating summary overlay...", file=sys.stderr)
    fig, ax = plt.subplots(figsize=(10, 7))

    for d in DIRS:
        log_dir = os.path.join(BASE_DIR, d, 'LOGS')
        if not os.path.isdir(log_dir):
            continue
        negs = scan_neg_profiles(log_dir)
        if not negs:
            continue
        for profile_num, model_number, neg_zones, fpath in negs:
            prof = mr.MesaData(fpath)
            logT = np.array(prof.logT, dtype=float)
            logP = np.array(prof.logP, dtype=float)
            gradT = np.array(prof.gradT, dtype=float)
            mask = gradT < 0
            ax.plot(logP[mask], logT[mask], '.', markersize=2, alpha=0.5,
                   label=utils.parse_name_for_plots(d) if d not in [x[0] for x in ax.get_legend_handles_labels()[1]] else "")

    ax.set_xlabel(utils.pretty_axis_name('logP'))
    ax.set_ylabel(utils.pretty_axis_name('logT'))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=6, loc='best', markerscale=3)
    ax.set_title('All negative gradT regions (logT vs logP)')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'summary_all_neg_gradT.png'), dpi=300)
    plt.close()
    print(f"Plots saved to {out_dir}/", file=sys.stderr)


if __name__ == '__main__':
    plot_dir()
