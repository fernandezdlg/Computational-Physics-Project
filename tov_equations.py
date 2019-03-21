# polytopic_eos.py contains function that returns the value of the specific internal energy of a polytopic fluid
# outlined as e = P / (y - 1)*rho , where P is some constant, gamma is the polytopic constant

from scipy import *
from polytropic_eos import *
from tov_int_methods import *
from convert_SI_to_cGM import *
import matplotlib.pyplot as plt
 
# These 'dpdr', 'dmdr' and 'dphidr' are in SI units
"""
def dpdr(rho,epsilon,r,p,m):
    return -G*(rho*(1+epsilon/c**2)+p/c**2)*(m+4*pi*r**3*p/c**2)/(r*(r-(2*G*m/c**2)))
    
def dmdr(r,rho,epsilon):
    return 4*pi*r**2*rho*(1+epsilon/(c**2))
    
def dphidr(m,r,p):
    return (m+4*pi*r**3*p/c**2)/(r*(r-2*G*m/c**2))
"""

# These 'dpdr', 'dmdr', 'dphidr' and 'dmbdr' are in cGM units
def dpdr(rho,epsilon,r,p,m):
    return -(rho*(1+epsilon)+p)*(m + 4*pi*(r**3)*p)/(r*(r-(2*m)))
   
def dmdr(r,rho,epsilon):
    return 4*pi*(r**2)*rho*(1+epsilon)

def dphidr(r,m,p):
    return -(m+(4*pi*(r**3)*p))/(r*(r-(2*m)))
    
def dmBdr(r,rho,m): 
    return 4*pi*(r**2)*rho / sqrt(1 - 2*m/r)
    
def main():
    # Rmax is set in meters
    # convert to cGM units
    Rmax = 15000
    Rmax = convert_SI_Length(Rmax)

    N = 100000
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
    
    # initial pressure value is in cGM units
    # Likewise mass and potential also
    pressure[0] = pressure_eos(rho_c, K, gamma)
    mass[0], potential[0] = 0., 0.
    
    
    
    ###############################################################
    ### Work in progress ###
    # Attempting to add Euler integration function to separate file to clean up code
    #pressure, mass = euler_int(radii, pressure, mass, dr, K, gamma)
    ###############################################################
    
    
    
    # Euler method  
    """
    for n,r_n in enumerate(radii):
        if n != 0:
            p_n = pressure[n-1]
            m_n = mass[n-1]
            phi_n = potential[n-1]
        
            rho_n = rho_eos(p_n, K, gamma)
            epsilon_n = sp_energy_eos(p_n, rho_n, gamma) 
        
            pressure[n] = p_n + (dr * dpdr(rho_n, epsilon_n, r_n, p_n, m_n))
            mass[n] = m_n + (dr * dmdr(r_n, rho_n, epsilon_n))
            potential[n] = phi_n + (dr * dphidr(r_n, m_n, p_n))
    """
    
    # RK2 method
    """
    for n,r_n in enumerate(radii):
        if n != 0:
            p_n = pressure[n-1]
            m_n = mass[n-1]
            phi_n = potential[n-1]
            
            rho_n = rho_eos(p_n, K, gamma)
            epsilon_n = sp_energy_eos(p_n, rho_n, gamma)
        
            p_k1 = dpdr(rho_n, epsilon_n, r_n, p_n, m_n)
            p_k2 = dpdr(rho_n, epsilon_n, r_n + 0.5*dr, p_n + 0.5*dr*p_k1, m_n)
            pressure[n] = p_n + dr*p_k2
    
            m_k1 = dmdr(r_n, rho_n, epsilon_n)
            m_k2 = dmdr(r_n + 0.5*dr, rho_n, epsilon_n)
            mass[n] = m_n + dr*m_k2
            
            phi_k1 = dphidr(r_n, m_n, p_n)
            phi_k2 = dphidr(r_n + 0.5*dr, m_n, p_n)
            potential[n] = phi_n + dr*phi_k2
    """
    
    # RK3 method
    """
    for n,r_n in enumerate(radii):
        if n != 0:
            p_n = pressure[n-1]
            m_n = mass[n-1]
            phi_n = potential[n-1]
            
            rho_n = rho_eos(p_n, K, gamma)
            epsilon_n = sp_energy_eos(p_n, rho_n, gamma)
        
            p_k1 = dpdr(rho_n, epsilon_n, r_n, p_n, m_n)
            p_k2 = dpdr(rho_n, epsilon_n, r_n + 0.5*dr, p_n + 0.5*dr*p_k1, m_n)
            p_k3 = dpdr(rho_n, epsilon_n, r_n + dr, p_n - dr*p_k1 + 3*dr*p_k2, m_n)
            pressure[n] = p_n + (dr/6.)*(p_k1 + 4*p_k2 + p_k3)
    
            m_k1 = dmdr(r_n, rho_n, epsilon_n)
            m_k2 = dmdr(r_n + 0.5*dr, rho_n, epsilon_n) 
            m_k3 = dmdr(r_n + dr, rho_n, epsilon_n)
            mass[n] = m_n + (dr/6.)*(m_k1 + 4*m_k2 + m_k3)
            
            phi_k1 = dphidr(r_n, m_n, p_n)
            phi_k2 = dphidr(r_n + 0.5*dr, m_n, p_n)
            phi_k3 = dphidr(r_n + dr, m_n, p_n)
            potential[n] = phi_n + (dr/6.)*(phi_k1 + 4*phi_k2 + phi_k3)
    """       
    
    # RK4 method
    
    for n,r_n in enumerate(radii):
        if n != 0:
            p_n = pressure[n-1]
            m_n = mass[n-1]
            phi_n = potential[n-1]
            
            rho_n = rho_eos(p_n, K, gamma)
            epsilon_n = sp_energy_eos(p_n, rho_n, gamma)
        
            p_k1 = dpdr(rho_n, epsilon_n, r_n, p_n, m_n)
            p_k2 = dpdr(rho_n, epsilon_n, r_n + 0.5*dr, p_n + 0.5*dr*p_k1, m_n)
            p_k3 = dpdr(rho_n, epsilon_n, r_n + 0.5*dr, p_n + 0.5*dr*p_k2, m_n)
            p_k4 = dpdr(rho_n, epsilon_n, r_n + dr, p_n + dr*p_k3, m_n)
            pressure[n] = p_n + (dr/6.)*(p_k1 + 2*p_k2 + 2*p_k3 + p_k4)
    
            m_k1 = dmdr(r_n, rho_n, epsilon_n)
            m_k2 = dmdr(r_n + 0.5*dr, rho_n, epsilon_n)
            m_k3 = dmdr(r_n + 0.5*dr, rho_n, epsilon_n)
            m_k4 = dmdr(r_n + dr, rho_n, epsilon_n)
            mass[n] = m_n + (dr/6.)*(m_k1 + 2*m_k2 + 2*m_k3 + m_k4)
            
            ### Not sure if you have to put p_n + 0.5**p_k1 terms???
            phi_k1 = dphidr(r_n, m_n, p_n)
            phi_k2 = dphidr(r_n + 0.5*dr, m_n, p_n)
            phi_k3 = dphidr(r_n + 0.5*dr, m_n, p_n)
            phi_k4 = dphidr(r_n + dr, m_n, p_n)
            potential[n] = phi_n + (dr/6.)*(phi_k1 + 2*phi_k2 + 2*phi_k3 + phi_k4)
    
   
    surfaceindex = where(pressure < 0)[0][0]
    radiusOfStar = (1 / convert_SI_Length(1))*radii[surfaceindex]
    massOfStar = (1 / convert_SI_Mass(1))*mass[surfaceindex]
    print('The radius of the star (meters): '+str(radiusOfStar))
    print('The mass of the star (kg): '+str(massOfStar))
    
    ### Unsure about the units of the potential and possibly comparing values in CGM units to 
    ### Values in SI??
    """
    phiOfStar = potential[surfaceindex]
    print('The potential of the star (?): '+str(phiOfStar))
    schwarz_metric = 0.5 * log(1 - (2*mass[surfaceindex] / radii[surfaceindex]))
    print('The Schwarzschild metric potential of the star (?): '+str(schwarz_metric))
    """
    
    massBaryon = zeros(N)
    
    for n, r_n in enumerate(radii):
        if n != 0:
            mB_n = massBaryon[n-1]
            m_n = mass[n-1]
            p_n = pressure[n-1]
            rho_n = rho_eos(p_n, K, gamma)
        
            mB_k1 = dmBdr(r_n, rho_n, massOfStar)
            mB_k2 = dmBdr(r_n + 0.5*dr, rho_n, m_n)
            mB_k3 = dmBdr(r_n + 0.5*dr, rho_n, m_n)
            mB_k4 = dmBdr(r_n + dr, rho_n, m_n)
        
            massBaryon[n] = mB_n + (dr/6.)*(mB_k1 + 2*mB_k2 + 2*mB_k3 + mB_k4)
        
    baryonicMassOfStar = (1 / convert_SI_Mass(1)) * massBaryon[surfaceindex]
    print('The baryonic mass of the star (kg): '+str(baryonicMassOfStar))
    
    plt.close()
     
    plt.subplots(1,3)
    plt.subplot(131)
    plt.plot(radii,pressure)
    plt.xlabel('radius (cGM units)'), plt.ylabel('pressure (cGM units)')
    
    plt.subplot(132)
    plt.plot(radii,mass)
    plt.plot(radii, massBaryon)
    plt.xlabel('radius (cGM units)'), plt.ylabel('mass (cGM units)')
    
    plt.subplot(133)
    plt.plot(radii,potential)
    plt.xlabel('radius (cGM units)'), plt.ylabel('potential (cGM units)')
    
    plt.show()
        
        
    
main()
        
        
        