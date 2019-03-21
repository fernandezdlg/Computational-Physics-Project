from polytropic_eos import * 
from tov_equations import *

def euler_int(radii, pressure, mass, dr, K, gamma):
    for n,r_n in enumerate(radii):
        if n != 0:
            p_n = pressure[n-1]
            m_n = mass[n-1]
        
            rho_n = rho_eos(p_n, K, gamma)
            epsilon_n = sp_energy_eos(p_n, rho_n, gamma) 
        
            pressure[n] = p_n + (dr * dpdr(rho_n, epsilon_n, r_n, p_n, m_n))
            mass[n] = m_n + (dr * dmdr(r_n, rho_n, epsilon_n))
            
    return pressure, mass
    
    
    
    
    