import numpy as np
from src.constants import pi
from src.tov_equations import tov_dMdr, tov_dPdr
from src.integrators import RK4

def TOV_solver(r0, P_c, h, eos, abs_tol=1e-6, rel_tol=1e-10):
    '''
    Implementation of a TOV solver

    -Takes:
        r0: initial radius in meters
        P_c: central pressure in SI units
        h: step size in meters
        eos: EoS Table
        abs_tol: absolute pressure tolerance in SI units
        rel_tol: relative pressure tolerance in SI units
    -Returns: Neutron star radius and mass for given central pressure
    '''
    # -- Initialize variables -- #
    eps_c = eos.eps_of_P(P_c)
    M_c = (4/3) * pi * (r0**3) * eps_c
    #P_stop = max(abs_tol, rel_tol * P_c)
    P_stop = rel_tol * P_c
    #P_stop = abs_tol

    # -- Compute first iteration with initial conditions -- #
    r, P, M = RK4(tov_dPdr, tov_dMdr, h, eos, r0, P_c, M_c)
    first_solution = (r, P, M)
    solutions = [first_solution]

    # -- Compute subsequent iterations until surface -- #
    max_steps = 15000
    steps = 0
    while P > P_stop and steps < max_steps:
        r, P, M = RK4(tov_dPdr, tov_dMdr, h, eos, r, P, M)
        solution = (r, P, M)
        solutions.append(solution)
        steps += 1
    if any(val < 0 for val in solutions[-1]) or not all(np.isfinite(val) for val in solutions[-1]):
        return None
    else:
        # -- Linear interpolation routine to estimate actual values for mass and radius of the star -- #
        P_final = 0.00
        r2, P2, M2 = solutions[-1]
        r1, P1, M1 = solutions[-2]
        r_final = r1 + ((P_final - P1) / (P2 - P1)) * (r2 - r1)
        M_final = M1 + ((r_final - r1) / (r2 - r1)) * (M2 - M1)
        return r_final, M_final