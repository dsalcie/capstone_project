import os
from pathlib import Path

BASE = Path(__file__).parent # Project directory

import numpy as np
from src.constants import r0, h, M_solar, rho, MEVFM3_TO_SI_PRESSURE
from src.eos import load_eos_data, TabulatedEOS
from src.solver import TOV_solver
from src.sequence import generate_star_sequence
from src.plotting import plot_data

Pc_min = 4 # MeVfm3
Pc_max = 1250 # MeVfm3
Number_of_stars = 700

EOS_file = 'eos_apr_apr_unified_crust.table' # <- modify this line for desired EOS file
data_folder = BASE / 'data'
file_path = data_folder / EOS_file

eps_SI, P_SI = load_eos_data(file_path)
eos = TabulatedEOS(eps_SI, P_SI)

def tabulate_star_sequence(EOS_file, r0, h, Pc_min, Pc_max, num, solver, eos):    
    name, ext = os.path.splitext(EOS_file)
    output_dir = BASE / 'outputs'
    output_file = output_dir / f'star_sequence_{name}.txt'

    star_sequence = generate_star_sequence(
        r0, h, 
        Pc_min, 
        Pc_max, 
        num, 
        solver,
        eos
        )

    # Generate tables of values for star sequence
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
        Pc_MeVfm3_range = np.linspace(Pc_min, Pc_max, num)
        for star, Pc in zip(star_sequence, Pc_MeVfm3_range):
            if star is not None:
                r, P, M, eps = star
                R_km = r / 1000 # m to km
                Pc_SI = Pc * MEVFM3_TO_SI_PRESSURE # convert to SI for eps_c
                M_M_solar = M / M_solar
                eps_c_SI = eos.cs_eps_of_P(Pc_SI) # interpolation in SI
                eps_c_MeVfm3 = eps_c_SI / rho # convert back to MeVfm3
                rows.append((Pc, R_km, M_M_solar, eps_c_MeVfm3))

        # Sort by radius (km) before writing.
        rows.sort(key=lambda row: row[2])
        for Pc_MeVfm3, R_km, M_M_solar, eps_c in rows:
            row = (
                f'{Pc_MeVfm3:15.8e}  {R_km:15.8e}   {M_M_solar:15.8e}   {eps_c:15.8e}\n'
            )
            file.write(row)
    print(f'Solutions have been printed in {output_file}')
    return output_file

if __name__ == "__main__":
    tabulate_star_sequence(EOS_file, r0, h, Pc_min, Pc_max, Number_of_stars, TOV_solver, eos)
    plot_data(BASE)
