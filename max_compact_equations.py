from scipy import *
from polytropic_eos import *
from convert_SI_to_cGM import *
import matplotlib.pyplot as plt

### Maximally compact EOS ###
# cGM units throughout
def dqdx(x,y,q):
    return -(y + 4*pi*q*(x**3))*(1+(2*q)) / (x*(x - (2*y)))

def dydx(x,q):
    return 4*pi*(x**2)*(1+q)
    
def compact_RK4(x_vals, q_vals, y_vals):
    dx = x_vals[1]-x_vals[0]
    
    for n,x_n in enumerate(x_vals):
        if n != 0:
            q_n = q_vals[n-1]
            y_n = y_vals[n-1]
            
            
            q_k1 = dqdx(x_n, y_n, q_n)
            q_k2 = dqdx(x_n + 0.5*dx, y_n, q_n + 0.5*dx*q_k1)
            q_k3 = dqdx(x_n + 0.5*dx, y_n, q_n + 0.5*dx*q_k2)
            q_k4 = dqdx(x_n + dx, y_n, q_n + dx*q_k3)
            q_vals[n] = q_n + (dx/6.)*(q_k1 + 2*q_k2 + 2*q_k3 + q_k4)
    
            y_k1 = dydx(x_n, q_n)
            y_k2 = dydx(x_n + 0.5*dx, q_n)
            y_k3 = dydx(x_n + 0.5*dx, q_n)
            y_k4 = dydx(x_n + dx, q_n)
            y_vals[n] = y_n + (dx/6.)*(y_k1 + 2*y_k2 + 2*y_k3 + y_k4)
    
    return q_vals, y_vals
    
def compact_main():
    # number of points
    N = 2000
    
    # conditions set in dimensionless units
    # in maximally compact configuration
    # q values are confined by q_max <= q <= 0
    q_max = 2.026
    
    q_vals, y_vals = zeros(N), zeros(N)
    q_vals[0] = q_max
    
    x_bound = 1.0
    x_vals = linspace(0, x_bound, N)
    
    
    
    # maximum e_0 value
    N_energy = 200
    e_0_max = 6 * 10**17
    e_0_vals = linspace(0, e_0_max, N_energy)
        
    maxmass_vals = zeros(N_energy)
    radius_vals = zeros(N_energy)
    
    for i, e_0 in enumerate(e_0_vals):
        #q_max = 2.026
    
        q_vals, y_vals = zeros(N), zeros(N) 
        q_vals[0] = q_max
    
        x_bound = 1
        x_vals = linspace(0, x_bound, N)
        
        q_vals, y_vals = compact_RK4(x_vals, q_vals, y_vals)
        
        # index of the point at the surface of the star
        surfaceindex = where(q_vals < 0)[0][0]
        
        # radius and mass in cGM units
        maxmass = y_vals[surfaceindex] * (c**4) / (G**3 * e_0)**0.5
        maxmass = (y_vals[surfaceindex] / (e_0*convert_SI_Density(1))**0.5)*(1/convert_SI_Mass(1))
        radius = (x_vals[surfaceindex] / (e_0*convert_SI_Density(1))**0.5)*(1/convert_SI_Length(1))
        
        maxmass_vals[i] = maxmass
        radius_vals[i] = radius
        
    chosenindex = where(radius_vals < 13800)[0][0]
    print('Maximum mass (with the same radius as the polytropic EoS with standard rho_c, K, gamma) is (kg): '+str(maxmass_vals[chosenindex]))    
        
    plt.close('all')
    
    #plt.subplots(1,2)
    #plt.subplot(121)
    plt.figure()
    plt.plot(e_0_vals, maxmass_vals)
    plt.xlabel('Surface energy density [J m$^{-2}$]'), plt.ylabel('Mass [kg]')
    plt.title('Total mass of neutron star')
    plt.tight_layout()
    plt.show()
    plt.savefig('maxMass.png', format='png', dpi=300)
    
    
    #plt.subplot(122)
    plt.figure()
    plt.plot(e_0_vals, radius_vals)
    plt.xlabel('Surface energy density [J m$^{-2}$]'), plt.ylabel('Radius [m]')
    plt.title('Radius of neutron star')
    plt.tight_layout()
    plt.show()
    plt.savefig('maxRadius.png', format='png', dpi=300)


compact_main()
    