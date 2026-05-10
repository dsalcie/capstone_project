import numpy as np
from scipy.interpolate import CubicSpline as cs
from src.constants import rho

def load_eos_data(file_path):   
    
    # Load data from EOS file
    data = np.loadtxt(file_path, skiprows=1)
    eps_MeVfm3 = data[:, 4]
    P_MeVfm3 = data[:, 3]

    # Unit conversion of data from MeV/fm^3 to SI units
    eps_SI = rho * eps_MeVfm3 # kg/m^3
    P_SI = P_MeVfm3 * 1.602176634e32 # J/m^3
    return eps_SI, P_SI

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
        self.cs_eps_of_P = cs(P_final, eps_final, bc_type='natural')

    def eps_of_P(self, P):
        P_clamped = np.clip(P, self.P_min, self.P_max)
        return self.cs_eps_of_P(P_clamped)