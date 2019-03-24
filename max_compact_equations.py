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
    
    # x values are confined from 0 (at centre of star) to minimum radius (at maximum compactness)
    x_min = 0.2404
    x_vals = linspace(0, x_min, N)
    
    q_vals, y_vals = zeros(N), zeros(N)
    q_vals[0] = q_max
    
    q_vals, y_vals = compact_RK4(x_vals, q_vals, y_vals)
    
    print(q_vals)
    print(y_vals)
    
    plt.close()
    
    plt.subplots(1,2)
    plt.subplot(121)
    plt.plot(x_vals,q_vals)
    plt.xlabel('x (dimensionless)'), plt.ylabel('q (dimensionless)')
    
    plt.subplot(122)
    plt.plot(x_vals,y_vals)
    plt.xlabel('x (dimensionless)'), plt.ylabel('y (dimensionless)')
    
    plt.show()
    

compact_main()
    