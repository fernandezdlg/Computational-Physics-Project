from scipy import *
from polytropic_eos import *
from convert_SI_to_cGM import *

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
    return (m+(4*pi*(r**3)*p))/(r*(r-(2*m))) #*-1
    
def dmBdr(r,rho,m): 
    return 4*pi*(r**2)*rho / sqrt(1 - 2*m/r)

### Newtonian equations ###
# cGM units used thoroughout
def dpdr_newton(rho,r,m):
    return -(rho*m)/(r**2)

def dmdr_newton(r, rho):
    return 4*pi*(r**2)*rho

def dphidr_newton(r,m,p):
    return m/(r**2)
