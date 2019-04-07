# polytopic_eos.py contains function that returns the value of the specific internal energy of a polytopic fluid
# outlined as e = P / (y - 1)*rho , where P is some constant, gamma is the polytopic constant

from scipy import *
from polytropic_eos import *
from dydt import *
from convert_SI_to_cGM import *
import matplotlib.pyplot as plt
 

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
            #mB_k3 = dmBdr(r_n + 0.5*dr, rho_n, massOfStar)
            mB_k4 = dmBdr(r_n + dr, rho_n, massOfStar)
        
            massBaryon[n] = mB_n + (dr/6.)*(mB_k1 + 4*mB_k2 +mB_k4)# + 2*mB_k3 + mB_k4)
    return massBaryon

def RK4_potential(radii, potential, mass, pressure, K, gamma):
    dr = radii[1]-radii[0]
    N = len(potential)
    
    for n in range(N-1, 0, -1):
        r_n = radii[n]
        p_n = pressure[n]
        m_n = mass[n]
        phi_n = potential[n]
            
        rho_n = rho_eos(p_n, K, gamma)
        epsilon_n = sp_energy_eos(p_n, rho_n, gamma)
    
        phi_k1 = dphidr(r_n, m_n, p_n)
        phi_k2 = dphidr(r_n - 0.5*dr, m_n, p_n)
        #phi_k3 = dphidr(r_n - 0.5*dr, m_n, p_n)
        phi_k4 = dphidr(r_n - dr, m_n, p_n)
    
        potential[n-1] = phi_n - (dr/6.)*(phi_k1 + 4*phi_k2 + phi_k4)# + 2*phi_k3
    
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
