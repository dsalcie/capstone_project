import os
from pathlib import Path

# -- Project directory -- #
BASE = Path(__file__).parent

import numpy as np
from src.constants import M_solar, rho, MEVFM3_TO_SI_PRESSURE
from src.eos import load_eos_data, TabulatedEOS
from src.solver import TOV_solver
from src.sequence import generate_star_sequence
from src.plotting import plot_data

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

# -- Two functions for generating tables -- #
def tabulate_star_sequence(EOS_file, r0, h, Pc_min, Pc_max, num, star_sequence, Pc_range, eos):    
    name, ext = os.path.splitext(EOS_file)
    output_dir = BASE / 'outputs'
    output_file = output_dir / f'star_sequence_{name}.txt'

    with open(output_file, 'w') as file:
        header = (
            'Stellar Sequence\n'
            f'Central Pressure Range = {Pc_min} to {Pc_max} MeVfm^3\n'
            f'Number of Stars = {num}\n'
            f'h = {h} meters\n'
            f'r0 = {r0} meters\n'
            f'{"P_MeV/fm3":>15}  '
            f'{"R_km":>15}   '
            f'{"M_M_solar":>15}    '
            f'{"Energy Density":>15}'
        )
        file.write(header + '\n')
        rows = []
        for star, Pc in zip(star_sequence, Pc_range):
            if star is not None:
                r, M, r2, M2 = star
                R_km = r2 / 1000 # m to km
                Pc_SI = Pc * MEVFM3_TO_SI_PRESSURE # convert to SI for eps_c
                M_M_solar = M2 / M_solar
                eps_c_SI = eos.eps_of_P(Pc_SI) # interpolation in SI
                eps_c_MeVfm3 = eps_c_SI / rho # convert back to MeVfm3
                rows.append((Pc, R_km, M_M_solar, eps_c_MeVfm3))

        for Pc_MeVfm3, R_km, M_M_solar, eps_c in rows:
            row = (
                f'{Pc_MeVfm3:15.8e}  {R_km:15.8e}   {M_M_solar:15.8e}   {eps_c:15.8e}\n'
            )
            file.write(row)
    print(f'Solutions have been printed in star_sequence_{name}')
    return output_file

def tov_solver_data_comparison(EOS_file, r0, h, Pc_min, Pc_max, num, star_sequence, Pc_range):
    '''
    Generate table to compare results from TOV Equation solver, where one column displays results that used 
    linear extrapolation routine and another column displays results that did not use linear extrapolation routine 
    to estimate the actual mass and radius of the neutron star
    '''
    name, ext = os.path.splitext(EOS_file)
    output_dir = BASE / 'outputs' / 'results comparison'
    output_file = output_dir / f'{name}.txt'
    
    with open(output_file, 'w') as file:
        header = (
            'Stellar Sequence\n'
            f'Central Pressure Range = {Pc_min} to {Pc_max} MeVfm^3\n'
            f'Number of Stars = {num}\n'
            f'h = {h} meters\n'
            f'r0 = {r0} meters\n'
            f'{"P_MeV/fm3":>15}        '
            f'{"R_km using LE":>15}     '
            f'{"R_km not using LE":>15} '
            f'{"M_M_solar using LE":>15} '
            f'{"M_M_solar not using LE":>15}'
        )
        file.write(header + '\n')
        rows = []
        for star, Pc in zip(star_sequence, Pc_range):
            if star is not None:
                r1, M1, r2, M2 = star
                R_km_LE = r1 / 1000 # m to km
                M_M_solar_LE = M1 / M_solar
                R_km_not_LE = r2 / 1000
                M_M_solar_not_LE = M2 / M_solar
                rows.append((Pc, R_km_LE, R_km_not_LE, M_M_solar_LE, M_M_solar_not_LE))
        for Pc_MeVfm3, R_km_LE, R_km_not_LE, M_M_solar_LE, M_M_solar_not_LE in rows:
            row = (
                f'{Pc_MeVfm3:15.8e}  {R_km_LE:15.8e}   {R_km_not_LE:15.8e}   {M_M_solar_LE:15.8e}   {M_M_solar_not_LE:15.8e}\n'
            )
            file.write(row)
    
    print(f'Results have been printed in {name} for comparison.')
    return output_file


if __name__ == "__main__":
    star_sequence = generate_star_sequence(
        r0, h, Pc_min, Pc_max, Number_of_stars, TOV_solver, eos
        )
    tabulate_star_sequence(
        EOS_file, r0, h, Pc_min, Pc_max, Number_of_stars, star_sequence, Pc_MeVfm3_range, eos
        )
    tov_solver_data_comparison(
        EOS_file, r0, h, Pc_min, Pc_max, Number_of_stars, star_sequence, Pc_MeVfm3_range
        )
    plot_data(BASE)
