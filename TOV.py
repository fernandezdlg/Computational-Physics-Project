# polytopic_eos.py contains function that returns the value of the specific internal energy of a polytopic fluid
# outlined as e = P / (y - 1)*rho , where P is some constant, y is the polytopic 

from polytropic_eos import *


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
    
