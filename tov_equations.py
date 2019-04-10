# polytopic_eos.py contains function that returns the value of the specific internal energy of a polytopic fluid
# outlined as e = P / (y - 1)*rho , where P is some constant, gamma is the polytopic constant
# int_methods contains the integrations methods used in this script

from scipy import *
from polytropic_eos import *
from int_methods import *
from dydt import *
from convert_SI_to_cGM import *
import matplotlib.pyplot as plt

#==============================================================================
# For Plotting
#==============================================================================
import matplotlib.pylab as pylab
# To plot with Serif font
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title(r'Neutron Stars Properties')
# Plotting parameters
params = {'legend.fontsize':'xx-large',
          'figure.figsize': (12, 6),
          'axes.labelsize': 'xx-large',
          'axes.titlesize': 'xx-large',
          'xtick.labelsize':'xx-large',
          'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)



def maximum_mass_calc(radii, N, K, gamma):
    # dr = difference in radius
    dr = radii[1]-radii[0]
    
    # maximum density of the star in kg m^-3
    N_rho = 1000
    rho_c_max = 101 * 10**17
    rho_c_min = 1 * 10**17
    rho_c_vals = linspace(rho_c_min, rho_c_max, N_rho)
    
    # convert to cGM units
    rho_c_vals = convert_SI_Density(1)*rho_c_vals
    
    pressure, mass = zeros(N), zeros(N)
    starMasses = zeros(N_rho)
    
    for i, rho_c in enumerate(rho_c_vals):
        pressure[0] = pressure_eos(rho_c, K, gamma)
        mass[0] = 0.0
        
        pressure, mass = RK4(radii, pressure, mass, K, gamma)
        
        # index of the point at the surface of the star
        # use index to find the mass of star
        surfaceindex = where(pressure < 0.)[0][0]
        massOfStar_cGM = mass[surfaceindex]
        
        starMasses[i] = massOfStar_cGM
        
    return rho_c_vals, starMasses


def poly_main():
    # Rmax is set in meters
    # convert to cGM units
    Rmax = 50000
    Rmax = convert_SI_Length(Rmax)
    
    N = 2000
    radii = linspace(0,Rmax,N)
    
    # culmulative sums of pressure, mass and potential
    # set as initial conditions which are sensible
    # rho_c is density at centre in kg m^-3
    rho_c = 5.0 * 10**17
    rho_c = convert_SI_Density(rho_c)
    
    # K and gamma are in cGM units already [This is independent of m or cm etc.]
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
    print('--------------------')
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
    
    # find optimal density for greatest max
    max_index = argmax(starMasses)
    opt_rho = rho_c_vals[max_index] * (1/convert_SI_Density(1))
    print('The maximum possible mass for K = '+str(K)+' and gamma = '+str(gamma)+' is (kg): '+str(max_mass))
    print('The maximum possible starting density (with same K, gamma parameters) is (kg m^-3): '+str(opt_rho))
    
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
    print('--------------------')
    
    plt.close()
    
    """
    plt.subplots(1,4)
    plt.subplot(141)
    plt.plot(radii,pressure, label=r'Relativistic (Polytropic) ($\gamma$='+str(gamma)+',K='+str(K)+')')
    plt.plot(radii,pressure_newton, label=r'Newton (Polytropic) ($\gamma$='+str(gamma)+',K='+str(K)+')')
    plt.xlabel('radius (cGM units)'),plt.ylabel('pressure (cGM units)'), plt.title('Pressure of neutron star')
    plt.legend()
    
    plt.subplot(142)
    plt.plot(radii,mass, label=r'Relativistic (Polytropic) ($\gamma$='+str(gamma)+',K='+str(K)+')')
    plt.plot(radii,massBaryon, label=r'Baryonic (Polytropic) ($\gamma$='+str(gamma)+',K='+str(K)+')')
    plt.plot(radii,mass_newton, label=r'Newton (Polytropic) ($\gamma$='+str(gamma)+',K='+str(K)+')')
    plt.xlabel('radius (cGM units)'), plt.ylabel('mass (cGM units)'), plt.title('Mass of neutron star')
    plt.legend()
    
    plt.subplot(143)
    plt.plot(radii[:surfaceindex],potential)
    plt.xlabel('radius (cGM units)'), plt.ylabel('potential (cGM units)')
    plt.title('Gravitational potential of neutron star ($\gamma$='+str(gamma)+',K='+str(K)+')')
    """
    
    #plt.subplot(144)
    plt.plot(rho_c_vals,starMasses)
    plt.xlabel(r'$\rho_c$ (cGM units)'), plt.ylabel('star mass (cGM units)')
    plt.title(r'Star masses for different densities ($\gamma$='+str(gamma)+',K='+str(K)+')')
    plt.tight_layout()
    plt.show()
    
          
poly_main()
    
    
    