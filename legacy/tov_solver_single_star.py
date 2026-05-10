import numpy as np
from scipy.interpolate import CubicSpline

# Fundamental Constants

c = 2.99792458e8 # m/s
G = 6.67430e-11 # m^3 kg^-1 s^-2
M_solar = 1.98847e30 # kg

# TOV Equations

def tov_dPdr(r, P, M):
    '''
    TOV equation: dP/dr
    '''
    eps = eos.cs_eps_of_P(P)
    return -(G * (eps + (P/c**2)) * (M + (4*np.pi*r**3 * (P/c**2))))/(r**2 * (1 - (2*G*M/(r*c**2))))

def tov_dMdr(r, P, M):
    '''
    TOV equation: dM/dr
    '''
    eps = eos.cs_eps_of_P(P)
    return 4 * (np.pi * r**2) * eps

# User input

EOS_file = input('Input EOS file: ')
P_c = float(input('Input central pressure: '))
Pc_SI = P_c * 1.602176634e32 # MeVfm3 to J/m^3

# Step 1: Load data from EOS.dat file

# Note: Ensure that EOS.dat files and tov_solver.py are in the same repository

data = np.loadtxt(EOS_file, max_rows=10000)
eps_MeVfm3 = data[:, 0]
P_MeVfm3 = data[:, 1]

# Step 2: Perform unit conversion from MeV/fm^3 to SI units
     
rho = 1.78266e15 # MeV/fm^3 to kg/m^3
eps_SI = rho * eps_MeVfm3 # kg/m^3
P_SI = P_MeVfm3 * 1.602176634e32 # J/m^3

# Step 3: Interpolate data in EOS file: ensure values are sorted and strictly increasing for CubicSpline.

sorted_eps_SI = np.sort(eps_SI)
sorted_P_SI = np.sort(P_SI)
mask_eps_SI = [sorted_eps_SI[i] # remove repeated values
               for i in range(len(sorted_eps_SI) - 1) 
               if sorted_eps_SI[i+1] - sorted_eps_SI[i] > 0]
mask_P_SI = [sorted_P_SI[i] # remove repeated values
               for i in range(len(sorted_P_SI) - 1) 
               if sorted_P_SI[i+1] - sorted_P_SI[i] > 0]

class TabulatedEOS:
    
    def __init__(self, eps_SI, P_SI):
        self.P_min, self.P_max = np.min(P_SI), np.max(P_SI)
        self.eps_min, self.eps_max = np.min(eps_SI), np.max(eps_SI)
        
        # Build splines 
        self.cs_eps_of_P = CubicSpline(P_SI, eps_SI, bc_type='natural')
        self.cs_P_of_eps = CubicSpline(eps_SI, P_SI, bc_type='natural')

    def eps_of_P(self, P):
        P_clamped = np.clip(P, self.P_min, self.P_max)
        return self.cs_eps_of_P(P_clamped)

    def P_of_eps(self, eps):
        eps_clamped = np.clip(eps, self.eps_min, self.eps_max)
        return self.cs_P_of_eps(eps_clamped)

eos = TabulatedEOS(mask_eps_SI, mask_P_SI)

# Step 4: Implement TOV Equation Solver to determine r, P(r), and M(r)

def RK4(f, g, h, t0, y0, z0):
    '''
    Implementation of Runge Kutta 4th Order-2nd ODE Solver, 1 iteration
    
    -Note: the TOV_solver function will take care of iterating through different r values, 
    hence Runge Kutta 4th Order will take care of integrating for the given t
    
    -Takes:
        t: independent variable
        f: f(t, y, z) 
        g: g(t, y, z)
        h: step
        t0: initial value of independent variable
        y0: y(0), initial condition
        z0: z(0), initial condition
    -Returns: t, y, z values
    '''
    # Initialize variables
    t_n = t0
    y_n = y0
    z_n = z0

    k_1y = f(t_n, y_n, z_n)
    k_1z = g(t_n, y_n, z_n)
    k_2y = f(t_n + h/2, y_n + ((h/2)*k_1y), z_n + ((h/2)*k_1z))
    k_2z = g(t_n + h/2, y_n + ((h/2)*k_1y), z_n + ((h/2)*k_1z))
    k_3y = f(t_n + h/2, y_n + ((h/2)*k_2y), z_n + ((h/2)*k_2z))
    k_3z = g(t_n + h/2, y_n + ((h/2)*k_2y), z_n + ((h/2)*k_2z))
    k_4y = f(t_n + h, y_n + (h * k_3y), z_n + (h * k_3z))
    k_4z = g(t_n + h, y_n + (h * k_3y), z_n + (h * k_3z))

    t_n += h
    y_n += (h/6)*(k_1y + (2*k_2y) + (2*k_3y) + k_4y)
    z_n += (h/6)*(k_1z + (2*k_2z) + (2*k_3z) + k_4z)
    return t_n, y_n, z_n

# Step 4: Solve TOV Equations using Runge-Kutta 4th Order ODE Solver to generate star profile

def TOV_solver(r0, P_c, h, abs_tol=1e-6, rel_tol=1e-8):
    '''
    Simple Implementation of a TOV solver

    -Takes:
        r0: initial radius (e.g., 10-100 m)
        h: step
        abs_tol: absolute tolerance
        rel_tol: relative tolerance
    -Returns: list of tuples containing r, P, M, and eps (Energy Density) values for chosen central pressure
    '''
    # Initialize variables
    eps_c = eos.cs_eps_of_P(P_c)
    M_c = (4/3) * np.pi * (r0**3) * eps_c
    P_stop = max(abs_tol, rel_tol * P_c)

    # Compute first iteration with initial conditions
    first_solution = RK4(tov_dPdr, tov_dMdr, h, r0, P_c, M_c) + (eps_c,)
    solutions = [first_solution]
    r, P, M, eps_c = solutions[0]

    # Compute subsequent iterations until surface
    iterations = 0
    while P > P_stop and iterations < 5000:
        r, P, M = RK4(tov_dPdr, tov_dMdr, h, r, P, M)
        eps = eos.eps_of_P(P)
        solution = (r, P , M, eps)
        if P < 0 or M < 0:
            iterations += 1
        else:
            solutions.append(solution)
            iterations += 1
    return solutions

r0 = 10 # in meters (m)
h = 5 # m

# Generate star_profile.txt file

output_file = f'star_profile_Pc_{P_c}.txt'
star_profile = TOV_solver(r0, Pc_SI, h)

# Print data from TOV solver
with open(output_file, 'w') as file:
    header = (
        f'{"R_km":>15}   '
        f'{"P_MeV/fm3":>15}   '
        f'{"M_M_solar":>15}   '
        f'{"Energy Density (ε)":>15}'
    )
    file.write(header + '\n')
    rows = []
    for i in range(len(star_profile) - 1):
        r, P, M, eps = star_profile[i]
        R_km = r / 1000 # m to km
        P_SI = P / 1.602176634e32 # convert back to MeVfm3
        M_M_solar = M / M_solar
        rows.append((R_km, P_SI, M_M_solar, eps))
    
    # Sort by radius (km) before writing.
    rows.sort(key=lambda row: row[0])
    for R_km, P_MeVfm3, M_M_solar, eps in rows:
        row = (
            f'{R_km:15.8e}   {P_MeVfm3:15.8e}   {M_M_solar:15.8e}   {eps:15.8e}\n'
        )
        file.write(row)

Mass_of_star = star_profile[-1][2] / M_solar
Radius_of_star = star_profile[-1][0] / 1000

print(f'Solutions have been printed in {output_file}')
print(f'Mass of Star (M_M_solar): {Mass_of_star}')
print(f'Radius of star (km): {Radius_of_star}')
