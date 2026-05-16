from src.constants import G, c, pi

def tov_dPdr(r, P, M, eos):
    '''
    TOV equation: dP/dr
    '''
    eps = eos.eps_of_P(P)
    return -(G * (eps + (P/c**2)) * (M + (4*pi*r**3 * (P/c**2))))/(r**2 * (1 - (2*G*M/(r*c**2))))

def tov_dMdr(r, P, M, eos):
    '''
    TOV equation: dM/dr
    '''
    eps = eos.eps_of_P(P)
    return 4 * (pi * r**2) * eps