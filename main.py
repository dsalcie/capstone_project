import os
from pathlib import Path

# -- Project directory -- #
BASE = Path(__file__).parent

import numpy as np
from src.eos import load_eos_data, TabulatedEOS
from src.solver import TOV_solver
from src.sequence import generate_star_sequence
from src.plotting import plot_data
from src.tabulate import tabulate_star_sequence, tov_solver_data_comparison

# -- Central pressure range -- #
Pc_min = np.float64(4.0) # MeVfm3
Pc_max = np.float64(1000.0) # MeVfm3
Number_of_stars = 1250

# -- Modify as needed for integrator routine -- #
h = 5 # m, step size
r0 = np.float64(10.0) # in meters, not zero to avoid singularity 

# -- Load EOS table -- #
EOS_file = 'eos_abht_qmc_rmf1_unified_crust.table' # <- modify this line for desired EOS file
data_folder = BASE / 'data'
file_path = data_folder / EOS_file

eps_SI, P_SI = load_eos_data(file_path)
eos = TabulatedEOS(eps_SI, P_SI)

# -- Compute Central Pressure range in MeVfm3 -- #
Pc_MeVfm3_range = np.linspace(Pc_min, Pc_max, Number_of_stars)

if __name__ == "__main__":
    star_sequence = generate_star_sequence(
        r0, h, Pc_min, Pc_max, Number_of_stars, TOV_solver, eos
        )
    tabulate_star_sequence(
        BASE, EOS_file, r0, h, Pc_min, Pc_max, Number_of_stars, star_sequence, Pc_MeVfm3_range, eos
        )
    tov_solver_data_comparison(
        BASE, EOS_file, r0, h, Pc_min, Pc_max, Number_of_stars, star_sequence, Pc_MeVfm3_range
        )
    plot_data(BASE)
