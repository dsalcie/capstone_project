import numpy as np
from scipy.interpolate import CubicSpline
import os
from pathlib import Path
import matplotlib.pyplot as plt

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

# Load data from EOS.dat file
EOS_file = 'eos_ds_cmf7.table' # <- modify this line for desired EOS file
BASE = Path(__file__).parent
data_folder = BASE / "data"
file_path = data_folder / EOS_file

data = np.loadtxt(file_path, skiprows=1)
eps_MeVfm3 = data[:, 4]
P_MeVfm3 = data[:, 3]

# Set parameters
Pc_min = 1.5 # MeVfm3
Pc_max = 425 # MeVfm3
Number_of_stars = 500

# Unit Conversion of Parameters
Pc_min_SI = Pc_min * 1.602176634e32 # MeVfm3 to J/m^3
Pc_max_SI = Pc_max * 1.602176634e32 # MeVfm3 to J/m^3

# Perform unit conversion of data from MeV/fm^3 to SI units  
rho = 1.78266e15 # MeV/fm^3 to kg/m^3
eps_SI = rho * eps_MeVfm3 # kg/m^3
P_SI = P_MeVfm3 * 1.602176634e32 # J/m^3

# Interpolate data in EOS file: ensure values are sorted and strictly increasing for CubicSpline.
class TabulatedEOS:
    
    def __init__(self, eps_SI, P_SI):
        self.P_min, self.P_max = np.min(P_SI), np.max(P_SI)
        self.eps_min, self.eps_max = np.min(eps_SI), np.max(eps_SI)
        
        # Sort by pressure
        idx = np.argsort(P_SI)
        P_sorted = P_SI[idx]
        eps_sorted = eps_SI[idx]

        # Remove duplicate pressures
        P_unique, unique_idx = np.unique(P_sorted, return_index=True)
        eps_unique = eps_sorted[unique_idx]

        # Enforce monotonic energy density
        mask = np.insert(np.diff(eps_unique) > 0, 0, True)

        P_final = P_unique[mask]
        eps_final = eps_unique[mask]

        # Build splines 
        self.cs_eps_of_P = CubicSpline(P_final, eps_final, bc_type='natural')

    def eps_of_P(self, P):
        P_clamped = np.clip(P, self.P_min, self.P_max)
        return self.cs_eps_of_P(P_clamped)

eos = TabulatedEOS(eps_SI, P_SI)

# Step 4: Implement TOV Equation Solver to determine r, P(r), and M(r)
def RK4(f, g, h, t0, y0, z0):
    '''
    Implementation of Runge Kutta 4th Order-2nd ODE Solver, 1 iteration
    
    -Note: the TOV_solver function will take care of iterating through different r values.
     Runge Kutta 4th Order will take care of integrating for the given t
    
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

# Solve TOV Equations using Runge-Kutta 4th Order ODE Solver to generate star profile
def TOV_solver(r0, P_c, h, abs_tol=1e-6, rel_tol=1e-8):
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
    M_c = (4/3) * np.pi * (r0**3) * eps_c
    P_stop = max(abs_tol, rel_tol * P_c)

    # Compute first iteration with initial conditions
    first_solution = RK4(tov_dPdr, tov_dMdr, h, r0, P_c, M_c) + (eps_c,)
    solutions = [first_solution]
    r, P, M, eps_c = solutions[0]

    # Compute subsequent iterations until surface
    while P > P_stop:
        r, P, M = RK4(tov_dPdr, tov_dMdr, h, r, P, M)
        eps = eos.eps_of_P(P)
        solution = (r, P , M, eps)
        solutions.append(solution)
    if any(val < 0 for val in solutions[-1]):
        return None
    else:
        return solutions[-1]

r0 = 5 # in meters (m)
h = 15 # m

# Generate star_sequence.txt file
name, ext = os.path.splitext(EOS_file)
output_dir = BASE / 'outputs'
output_file = output_dir / f'star_sequence_{name}.txt'

star_sequence = [TOV_solver(r0, P_c, h) 
                 for P_c in np.linspace(Pc_min_SI, Pc_max_SI, Number_of_stars)]

# Print data from TOV solver
with open(output_file, 'w') as file:
    header = (
        'Stellar Sequence\n'
        f'Central Pressure Range = {Pc_min} to {Pc_max} MeVfm^3\n'
        f'Number of Stars = {Number_of_stars}\n'
        f'h = {h} meters\n'
        f'r0 = {r0} meters\n'
        f'{"P_MeV/fm3":>15}  '
        f'{"R_km":>15}   '
        f'{"M_M_solar":>15}    '
        f'{"Energy Density":>15}'
    )
    file.write(header + '\n')
    rows = []
    Pc_MeVfm3_range = np.linspace(Pc_min, Pc_max, Number_of_stars)
    for i, Pc in zip(range(0, len(star_sequence) - 1), Pc_MeVfm3_range):
        if star_sequence[i] is not None:
            r, P, M, eps = star_sequence[i]
            R_km = r / 1000 # m to km
            Pc_SI = Pc * 1.602176634e32 # convert to SI for eps_c
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

plt.figure()
for out in output_dir.glob('star_sequence_eos*'):
    results = np.loadtxt(f'{out}', skiprows=6)
    R_vals = results[:,1]
    M_vals = results[:,2]
    
    # fix file name for legend
    stem = out.stem
    label_raw = stem.split('sequence_', 1)[-1]
    label_clean = label_raw.replace('_', ' ')
    label_pretty = label_clean.upper()
    
    plt.plot(R_vals, M_vals, marker=None, label=f'{label_pretty}')
plt.xlabel('R (km)')
plt.ylabel(r'M(M$_{\odot}$)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('all_eos.pdf', format='pdf', bbox_inches = 'tight', dpi=300)
plt.show()