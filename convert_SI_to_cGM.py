# file that contains all necessary functions to convert SI units to c = G = solar mass = 1 unit system
# ensure SI units used are kilograms for mass, metres for distance and seconds for time

# speed of light in vacuum as quoted by NIST CODATA in SI [m s^-1]
c = 2.99792458 * 10**8 
c_err = 0.0

# gravitational constant as quoted by NIST CODATA in SI [m^3 kg^-1 s^-2]
G = 6.67408 * 10**-11 
G_err = 0.00031 * 10**-11

# solar mass constant as quoted by NIST CODATA in SI [kg]
M0 = 1.98847 * 10**30  
M0_err = 0.00007 * 10**30

# mass is input in kg
def convert_SI_Mass(mass):
    return mass / M0

# time is input in seconds
def convert_SI_Time(time):
    return time * c**3 * G**-1 * M0**-1

# length is input in metres
def convert_SI_Length(length):
    return length * G**-1 * M0**-1 * c**2

# density is input in kg m^-3
def convert_SI_Density(density):
    return density * convert_SI_Length(1)**-3 * convert_SI_Mass(1)

# energy is input in kg m^2 s^-2
def convert_SI_Energy(energy):
    return energy * convert_SI_Mass(1) * convert_SI_Time(1)**-2 * convert_SI_Length(1)**2

# specific energy is input in m^2 s^-2
def convert_SI_SPEnergy(spenergy):
    return spenergy * convert_SI_Time(1)**-2 * convert_SI_Length(1)**2

# pressure is input in kg m^-1 s^-2
def convert_SI_Pressure(pressure):
    return pressure * convert_SI_Mass(1) * convert_SI_Length(1)**-1 * convert_SI_Time(1)**-2

# potential is input in m^2 s^-2
def convert_SI_Potential(potential):
    return potential * c**-2

"""
# test function used to check unit conversions are correct
def test():
    newTestVal = convert_SI_Pressure(1)
    print(newTestVal)
    
test()
"""