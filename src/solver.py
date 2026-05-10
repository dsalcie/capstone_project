from src.constants import pi
from src.tov_equations import tov_dMdr, tov_dPdr
from src.integrators import RK4

def TOV_solver(r0, P_c, h, eos, abs_tol=1e-6, rel_tol=1e-8):
    '''
    Simple Implementation of a TOV solver

    -Takes:
        r0: initial radius (e.g., 10-100 m)
        P_c: central pressure
        h: step
        abs_tol: absolute tolerance
        rel_tol: relative tolerance
    -Returns: list of tuples containing r, P, M, and eps (Energy Density) values for chosen central pressure
    '''
    # Initialize variables
    eps_c = eos.cs_eps_of_P(P_c)
    M_c = (4/3) * pi * (r0**3) * eps_c
    P_stop = max(abs_tol, rel_tol * P_c)

    # Compute first iteration with initial conditions
    first_solution = RK4(tov_dPdr, tov_dMdr, h, eos, r0, P_c, M_c) + (eps_c,)
    solutions = [first_solution]
    r, P, M, eps_c = solutions[0]

    # Compute subsequent iterations until surface
    while P > P_stop:
        r, P, M = RK4(tov_dPdr, tov_dMdr, h, eos, r, P, M)
        eps = eos.eps_of_P(P)
        solution = (r, P , M, eps)
        solutions.append(solution)
    if any(val < 0 for val in solutions[-1]):
        return None
    else:
        return solutions[-1]