# polytopic_eos.py contains function that returns the value of the specific internal energy of a polytopic fluid
# outlined as e = P / (y - 1)*rho , where P is some constant, y is the polytopic 

from scipy import *
from polytropic_eos import *
from int_methods import *
from convert_SI_to_cGM import *
import matplotlib.pyplot as plt

# Function that takes all quantities and converts them all to computational units 
"""
def dpdr(rho,epsilon,r,p,m):
    return -G*(rho*(1+epsilon/c**2)+p/c**2)*(m+4*pi*r**3*p/c**2)/(r*(r-(2*G*m/c**2)))
    
def dmdr(r,rho,epsilon):
    return 4*pi*r**2*rho*(1+epsilon/(c**2))
    
def dphidr(m,r,p):
    return (m+4*pi*r**3*p/c**2)/(r*(r-2*G*m/c**2))
"""

# These 'dpdr' and 'dmdr' are in cGM units
def dpdr(rho,epsilon,r,p,m):
    return -(rho*(1+epsilon)+p)*(m + 4*pi*(r**3)*p)/(r*(r-(2*m)))
   
def dmdr(r,rho,epsilon):
    return 4*pi*(r**2)*rho*(1+epsilon)
    
def Gmass(r,rho,epsilon): # integrand
    return 4*pi*r**2*rho*(1+epsilon/(c**2))
    
def Bmass(r,rho,M): # integrand, M: total mass
    #return 4*pi*r**2*rho*sqrt(det(g))
    return (sqrt(1-2*G*M/(r*c**2)))**(-1)*4*pi*r^2*rho
    
def main():
    # Rmax is set to 50km = 50,000m
    # convert to cGM units
    Rmax = 50000
    Rmax = convert_SI_Length(Rmax)

    N = 10000
    radii = linspace(0,Rmax,N)
    dr = radii[1]-radii[0]
    
    #ctr = [0,0,0]
    pressure, mass, potential = zeros(N), zeros(N), zeros(N)
    
    # culmulative sums of pressure, mass and potential
    # set as initial conditions which are sensible
    # rho_c is density at centre in kg m^-3
    rho_c = 5.0 * 10**17
    rho_c = convert_SI_Density(rho_c)
    
    # K and gamma are in cGM units already
    gamma = 2.75
    K = 30000
    
    """
    epsi0 = sp_energy_eos(p, K, gamma)
    print(epsi0)
    
    p_alt = pressure_eos(5.0e11, 1.98183e-6, gamma)
    epsi0_alt = sp_energy_eos(p_alt, 1.98183e-6, gamma)
    print(epsi0_alt)
    epsi0_alt_cGM = convert_SI_SPEnergy(epsi0_alt)
    print(epsi0_alt_cGM)
    """
    
    
    # initial pressure value is in cGM units
    # Likewise mass and potential also
    pressure[0] = pressure_eos(rho_c, K, gamma)
    mass[0], potential[0] = 0., 0.
    
    # Euler method  
    for i,r_n in enumerate(radii):
        if i != 0:
            p_n = pressure[i-1]
            m_n = mass[i-1]
        
            rho_n = rho_eos(p_n, K, gamma)
        
            epsilon_n = sp_energy_eos(p_n, rho_n, gamma) 
        
            pressure[i] = p_n + (dr * dpdr(rho_n, epsilon_n, r_n, p_n, m_n))
            mass[i] = m_n + (dr * dmdr(r_n, rho_n, epsilon_n))

    
    
    
    # RK4 method
    """
    for i, r in enumerate(radii[1:]):
        p_n, m_n = pressure[i], mass[i]
        rho_n = rho_eos(p_n, K, gamma)
        epsilon_n = sp_energy_eos(p_n, rho_n, gamma)
        
        k1 = dpdr(rho_n, epsilon_n, r, p_n, m_n)
        k2 = dpdr(rho_n, epsilon_n, r + 0.5*dr, p_n + 0.5*dr*k1, m_n)
        k3 = dpdr(rho_n, epsilon_n, r + 0.5*dr, p_n + 0.5*dr*k2, m_n)
        k4 = dpdr(rho_n, epsilon_n, r + dr, p_n + dr*k3, m_n)
        
        pressure[i+1] = p_n + (dr/6.)*(k1 + 2*k2 + 2*k3 + k4)
    """
    
   
    surfaceindex = where(pressure < 0)[0][0]
    radiusOfStar = (1 / convert_SI_Length(1))*radii[surfaceindex]
    massOfStar = (1 / convert_SI_Mass(1))*mass[surfaceindex]
    print('The radius of the star (meters): ' +str(radiusOfStar))
    print('The mass of the star (kg): ' +str(massOfStar))
    
    plt.close()
     
    plt.subplots(1,2)
    plt.subplot(121)
    plt.plot(radii,pressure)
    plt.xlabel('radius (cGM units)'), plt.ylabel('pressure (cGM units)')
    
    plt.subplot(122)
    plt.plot(radii,mass)
    plt.xlabel('radius (cGM units)'), plt.ylabel('mass (cGM units)')
    plt.show()
        
        
    
main()
        
        
        