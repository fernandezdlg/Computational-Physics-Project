# polytopic_eos.py contains function that returns the value of the specific internal energy of a polytopic fluid
# outlined as e = P / (y - 1)*rho , where P is some constant, gamma is the polytopic constant

from scipy import *
from polytropic_eos import *
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
### TOV equations ###
# cGM units used thoroughout
def dpdr(rho,epsilon,r,p,m):
    return -(rho*(1+epsilon)+p)*(m + 4*pi*(r**3)*p)/(r*(r-(2*m)))
   
def dmdr(r,rho,epsilon):
    return 4*pi*(r**2)*rho*(1+epsilon)

def dphidr(r,m,p):
    return -(m+(4*pi*(r**3)*p))/(r*(r-(2*m)))
    
def dmBdr(r,rho,m): 
    return 4*pi*(r**2)*rho / sqrt(1 - 2*m/r)

### Newtonian equations ###
# cGM units used thoroughout
def dpdr_newton(rho,r,m):
    return -(rho*m)/(r**2)

def dmdr_newton(r, rho):
    return 4*pi*(r**2)*rho

def dphidr_newton(r,m,p):
    return -m/(r**2)


### Integration methods ###
def euler(radii, pressure, mass, K, gamma):
    dr = radii[1]-radii[0]
    
    for n,r_n in enumerate(radii):
        if n != 0:
            p_n = pressure[n-1]
            m_n = mass[n-1]
        
            rho_n = rho_eos(p_n, K, gamma)
            epsilon_n = sp_energy_eos(p_n, rho_n, gamma) 
        
            pressure[n] = p_n + (dr * dpdr(rho_n, epsilon_n, r_n, p_n, m_n))
            mass[n] = m_n + (dr * dmdr(r_n, rho_n, epsilon_n))

    return pressure, mass

def RK2(radii, pressure, mass, K, gamma):
    dr = radii[1]-radii[0]
    
    for n,r_n in enumerate(radii):
        if n != 0:
            p_n = pressure[n-1]
            m_n = mass[n-1]
            
            rho_n = rho_eos(p_n, K, gamma)
            epsilon_n = sp_energy_eos(p_n, rho_n, gamma)
        
            p_k1 = dpdr(rho_n, epsilon_n, r_n, p_n, m_n)
            p_k2 = dpdr(rho_n, epsilon_n, r_n + 0.5*dr, p_n + 0.5*dr*p_k1, m_n)
            pressure[n] = p_n + dr*p_k2
    
            m_k1 = dmdr(r_n, rho_n, epsilon_n)
            m_k2 = dmdr(r_n + 0.5*dr, rho_n, epsilon_n)
            mass[n] = m_n + dr*m_k2
            
    return pressure, mass

def RK3(radii, pressure, mass, K, gamma):
    dr = radii[1]-radii[0]
    
    for n,r_n in enumerate(radii):
        if n != 0:
            p_n = pressure[n-1]
            m_n = mass[n-1]
            
            print(m_n)
            
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
            
    return pressure, mass

def RK4(radii, pressure, mass, K, gamma):
    dr = radii[1]-radii[0]
    
    for n,r_n in enumerate(radii):
        if n != 0:
            p_n = pressure[n-1]
            m_n = mass[n-1]
            
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
    
    return pressure, mass

def RK4_baryonic(radii, massBaryon, massOfStar, pressure, K, gamma):
    dr = radii[1]-radii[0]
    
    for n, r_n in enumerate(radii):
        if n != 0:
            mB_n = massBaryon[n-1]
            p_n = pressure[n-1]
            rho_n = rho_eos(p_n, K, gamma)
        
            mB_k1 = dmBdr(r_n, rho_n, massOfStar)
            mB_k2 = dmBdr(r_n + 0.5*dr, rho_n, massOfStar)
            mB_k3 = dmBdr(r_n + 0.5*dr, rho_n, massOfStar)
            mB_k4 = dmBdr(r_n + dr, rho_n, massOfStar)
        
            massBaryon[n] = mB_n + (dr/6.)*(mB_k1 + 2*mB_k2 + 2*mB_k3 + mB_k4)
            
    return massBaryon

def RK4_potential(radii, potential, mass, pressure, K, gamma):
    dr = radii[1]-radii[0]
    surfaceindex = len(potential)
    
    for n in range(surfaceindex-1, 0, -1):
        r_n = radii[n]
        p_n = pressure[n]
        m_n = mass[n]
        phi_n = potential[n]
            
        rho_n = rho_eos(p_n, K, gamma)
        epsilon_n = sp_energy_eos(p_n, rho_n, gamma)
    
        phi_k1 = dphidr(r_n, m_n, p_n)
        phi_k2 = dphidr(r_n - 0.5*dr, m_n, p_n)
        phi_k3 = dphidr(r_n - 0.5*dr, m_n, p_n)
        phi_k4 = dphidr(r_n - dr, m_n, p_n)
    
        potential[n-1] = phi_n - (dr/6.)*(phi_k1 + 2*phi_k2 + 2*phi_k3 + phi_k4)
    
    return potential

def Newton_RK4(radii, pressure, mass, K, gamma):
    dr = radii[1]-radii[0]
    
    for n,r_n in enumerate(radii):
        if n != 0:
            p_n = pressure[n-1]
            m_n = mass[n-1]
        
            rho_n = rho_eos(p_n, K, gamma)
        
            p_k1 = dpdr_newton(rho_n, r_n, m_n)
            p_k2 = dpdr_newton(rho_n, r_n + 0.5*dr, m_n)
            p_k3 = dpdr_newton(rho_n, r_n + 0.5*dr, m_n)
            p_k4 = dpdr_newton(rho_n, r_n + dr, m_n)
            pressure[n] = p_n + (dr/6.)*(p_k1 + 2*p_k2 + 2*p_k3 + p_k4)
            
            m_k1 = dmdr_newton(r_n, rho_n)
            m_k2 = dmdr_newton(r_n + 0.5*dr, rho_n)
            m_k3 = dmdr_newton(r_n + 0.5*dr, rho_n)
            m_k4 = dmdr_newton(r_n + dr, rho_n)
            mass[n] = m_n + (dr/6.)*(m_k1 + 2*m_k2 + 2*m_k3 + m_k4)
    
    return pressure, mass

def maximum_mass_calc(radii, N, K, gamma):
    # dr = difference in radius
    dr = radii[1]-radii[0]
    
    # maximum density of the star in kg m^-3
    N_rho = 100
    rho_c_max = 50 * 10**17
    rho_c_vals = linspace(rho_c_max/10, rho_c_max, N_rho)
    
    pressure, mass = zeros(N), zeros(N)
    starMasses = zeros(N_rho)
    
    for i, rho_c in enumerate(rho_c_vals):
        rho_c = convert_SI_Density(rho_c)
        
        pressure[0] = pressure_eos(rho_c, K, gamma)
        mass[0] = 0.0
        
        pressure, mass = RK4(radii, pressure, mass, K, gamma)
        
        # index of the point at the surface of the star
        # use index to find the mass of star
        surfaceindex = where(pressure < 0)[0][0]
        massOfStar_cGM = mass[surfaceindex]
        
        starMasses[i] = massOfStar_cGM
        
    return rho_c_vals, starMasses

    
def main():
    # Rmax is set in meters
    # convert to cGM units
    Rmax = 50000
    Rmax = convert_SI_Length(Rmax)

    N = 5000
    radii = linspace(0,Rmax,N)
    
    # culmulative sums of pressure, mass and potential
    # set as initial conditions which are sensible
    # rho_c is density at centre in kg m^-3
    rho_c = 5.0 * 10**17
    rho_c = convert_SI_Density(rho_c)
    
    # K and gamma are in cGM units already
    gamma = 2.75
    K = 30000 
    
    pressure, mass = zeros(N), zeros(N)
    
    # initial pressure value is in cGM units
    # Likewise mass and potential also
    pressure[0] = pressure_eos(rho_c, K, gamma)
    mass[0] = 0.  
    
    #pressure, mass = euler(radii, pressure, mass, K, gamma)
    #pressure, mass = RK2(radii, pressure, mass, K, gamma)
    #pressure, mass = RK3(radii, pressure, mass, K, gamma)
    pressure, mass = RK4(radii, pressure, mass, K, gamma)
    
    # index of the point at the surface of the star
    surfaceindex = where(pressure < 0)[0][0]
    
    # radius and mass in cGM units
    radiusOfStar_cGM = radii[surfaceindex]
    massOfStar_cGM = mass[surfaceindex]
    
    # radius and mass in SI units
    radiusOfStar = (1 / convert_SI_Length(1))*radiusOfStar_cGM
    massOfStar = (1 / convert_SI_Mass(1))*massOfStar_cGM
    print('The radius of the star (meters): '+str(radiusOfStar))
    print('The mass of the star (kg): '+str(massOfStar))
    
    
    potential = zeros(surfaceindex)
    potential[-1] = 0.5 * log(1 - (2*massOfStar_cGM / radiusOfStar_cGM))
    potential = RK4_potential(radii, potential, mass, pressure, K, gamma)

    nearCentrePotential = potential[1] * (1/convert_SI_Potential(1))
    print('The near centre of mass potential is (J/kg): '+str(nearCentrePotential))
    
    massBaryon = zeros(N)
    massBaryon = RK4_baryonic(radii, massBaryon, massOfStar_cGM, pressure, K, gamma)
    
    baryonicMassOfStar = (1/convert_SI_Mass(1)) * massBaryon[surfaceindex]
    print('The baryonic mass of the star (kg): '+str(baryonicMassOfStar))
    
    # find star masses for different starting densities
    rho_c_vals, starMasses = maximum_mass_calc(radii, N, K, gamma)
    max_mass = (1/convert_SI_Mass(1)) * max(starMasses)
    print('The maximum possible mass for K = '+str(K)+' and gamma = '+str(gamma)+' is (kg): '+str(max_mass))
    
    
    ##### Using Newtonian equations #####
    pressure_newton, mass_newton = zeros(N), zeros(N)
    pressure_newton[0] = pressure_eos(rho_c, K, gamma)
    mass_newton[0] = 0. 
    
    pressure_newton, mass_newton = Newton_RK4(radii, pressure_newton, mass_newton, K, gamma)
    
    # index of the point at the surface of the star
    surfaceindex_newton = where(pressure_newton < 0)[0][0]
    
    # radius and mass in cGM units
    radiusOfStar_newton_cGM = radii[surfaceindex_newton]
    massOfStar_newton_cGM = mass_newton[surfaceindex_newton]
    
    # radius and mass in SI units
    radiusOfStar_newton = (1 / convert_SI_Length(1))*radiusOfStar_newton_cGM
    massOfStar_newton = (1 / convert_SI_Mass(1))*massOfStar_newton_cGM
    print('The newton radius of the star (meters): '+str(radiusOfStar_newton))
    print('The newton mass of the star (kg): '+str(massOfStar_newton))
    
    plt.close()
    
    plt.subplots(1,4)
    plt.subplot(141)
    plt.plot(radii,pressure)
    plt.plot(radii,pressure_newton)
    plt.xlabel('radius (cGM units)'), plt.ylabel('pressure (cGM units)')
    
    plt.subplot(142)
    plt.plot(radii,mass)
    plt.plot(radii,massBaryon)
    plt.plot(radii,mass_newton)
    plt.xlabel('radius (cGM units)'), plt.ylabel('mass (cGM units)')
    
    plt.subplot(143)
    plt.plot(radii[:surfaceindex],potential)
    plt.xlabel('radius (cGM units)'), plt.ylabel('potential (cGM units)')
    
    plt.subplot(144)
    plt.plot(rho_c_vals,starMasses)
    plt.xlabel(r'$\rho_c$ (cGM units)'), plt.ylabel('star mass (cGM units)')
    
    plt.show()
        
        
    
main()
        
        
        