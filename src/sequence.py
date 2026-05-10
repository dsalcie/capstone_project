import numpy as np

def generate_star_sequence(r0, h, Pc_min_MeVfm3, Pc_max_MeVfm3, num, solver, eos):
    
    # Unit Conversion of Parameters
    Pc_min_SI = Pc_min_MeVfm3 * 1.602176634e32 # MeVfm3 to J/m^3
    Pc_max_SI = Pc_max_MeVfm3 * 1.602176634e32 # MeVfm3 to J/m^3
    
    return [
        solver(r0, P_c, h, eos) 
        for P_c in np.linspace(Pc_min_SI, Pc_max_SI, num)
        ]