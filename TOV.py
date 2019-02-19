# polytopic_eos.py contains function that returns the value of the specific internal energy of a polytopic fluid
# outlined as e = P / (y - 1)*rho , where P is some constant, y is the polytopic 

from scipy import *
from polytropic_eos import *
from int_methods import *
from convert_SI_to_cGM import *
import matplotlib.pyplot as plt

# Function that takes all quantities and converts them all to computational units
#def convert(G,c,epsilon,) 

def dpdr(rho,epsilon,p,m):
    return -G*(rho*(1+epsilon/c**2)+p/c**2)*(m+4*pi*r**3*p/c**2)/(r*(r-2*G*m/c**2))
    
def dmdr(r,rho,epsilon):
    return 4*pi*r**2*rho*(1+epsilon/(c**2))
    
def dphidr(m,r,p):
    return (m+4*pi*r**3*p/c**2)/(r*(r-2*G*m/c**2))
    
def Gmass(r,rho,epsilon): # integrand
    return 4*pi*r**2*rho*(1+epsilon/(c**2))
    
def Bmass(r,rho,M): # integrand, M: total mass
    #return 4*pi*r**2*rho*sqrt(det(g))
    return (sqrt(1-2*G*M/(r*c**2)))**(-1)*4*pi*r^2*rho
    
def main():
    Rmax = 50000
    N = 5001
    radii = linspace(0,Rmax,N)
    dr = Rmax/(N-1)
    
    ctr = [0,0,0]
    pressure, mass, potential = zeros(N), zeros(N), zeros(N)
    
    # culmulative sums of pressure, mass and potential
    # set as initial conditions which are sensible
    rho_c = 5.0 * 10**14
    gamma = 2.75
    K = 30000
    pressure[0] = pressure_eos(rho_c, K, gamma)
    mass[0], potential[0] = 0, 0
    
    for i,r in enumerate(radii[1:]):
        p_n = pressure[i]
        m_n = mass[i]
        rho_n = rho_eos(p_n, K, gamma)
        
        epsilon = sp_energy_eos(p_n, rho_n, gamma) 
        pressure[i+1] = p_n + (dr*dpdr(rho_n, epsilon, p_n, m_n))
    
    fig = plt.figure()    
    plt.plot(radii,pressure)
    plt.show()
        
        
    
main()
        
        
        