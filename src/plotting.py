import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def plot_data(BASE: Path):
    
    output_dir = BASE / 'outputs'
    plot_file = output_dir / 'all_eos.pdf'
    
    nicer_dir = BASE / 'data' / 'NICER_Posteriors'

    fig, ax = plt.subplots()

    # Plot NICER Posterior Distributions
    for nicer in nicer_dir.glob('PSR_J*.csv'):
        data = np.genfromtxt(nicer, delimiter=',')
        R_vals = data[:,0]
        M_vals = data[:,1]

        # fix file name for legend
        nicer_file_name = nicer.stem
        label_raw = nicer_file_name.replace('_', ' ')
        label_clean = label_raw.upper()

        ax.fill(
            R_vals,
            M_vals,
            label=f'{label_clean}',
            alpha=0.25,
            linewidth=1.2,
            zorder=1,
        )

    # Plot star sequences
    for out in output_dir.glob('star_sequence_eos*'):
        results = np.loadtxt(out, skiprows=6)
        R_vals = results[:,1]
        M_vals = results[:,2]

        # fix file name for legend
        output_file_name = out.stem
        label_raw = output_file_name.split('sequence_', 1)[-1]
        label_clean = label_raw.replace('_', ' ')
        label_pretty = label_clean.upper()

        ax.plot(R_vals, M_vals, marker=None, label=f'{label_pretty}', zorder=2)

    ax.set_xlabel('R (km)')
    ax.set_ylabel(r'M(M$_{\odot}$)')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(plot_file, format='pdf', bbox_inches='tight', dpi=300)
    plt.close(fig)
