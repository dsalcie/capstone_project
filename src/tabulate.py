import os
from constants import MEVFM3_TO_SI_PRESSURE, M_solar, rho

def tabulate_star_sequence(BASE, EOS_file, r0, h, Pc_min, Pc_max, num, star_sequence, Pc_range, eos):    
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

def tov_solver_data_comparison(BASE, EOS_file, r0, h, Pc_min, Pc_max, num, star_sequence, Pc_range):
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