import numpy as np
import matplotlib.pyplot as plt
#from math import Rational
from sympy import legendre, diff, integrate, factorial, exp, Rational,lambdify, Piecewise,simplify
from sympy.abc import x, t, u, z

bandwidth_data = np.loadtxt("bwoutput_scheme=euler_z0=0.50_N=100000000_h=2.0e-05.txt", delimiter="\t")

N_legendre = 50 # number of legendre polynomial terms
K = Piecewise( (Rational(3,4)*(1 - x**2), x*x<=1), (0, True) )
C = 1
for n in range(1,N_legendre): 
    C += (2*n+1) * legendre(n, 2*u-1) * legendre(n,0) * exp(-6*n*(n+1)*t)

bandwidth = 0.0007
ugrid = np.linspace(0,1.0, 51)
t_obs_list = [0.036, 0.072, 0.108, 0.144, 0.180, 0.216]
dt_list = [3e-3,1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]
r1 = np.zeros(( len(dt_list), len(t_obs_list),len(ugrid) ))
for i in range(len(dt_list)):
     for j in range(len(t_obs_list)):
        bandwidth = bandwidth_data[i,j]
        print(bandwidth)
        C_integrand = C * K.subs(x,(u-z)/bandwidth)
        for index in range(len(ugrid)):
            C_temp = C_integrand.subs({u: ugrid[index], t: t_obs_list[j]})
            r1[i, j, index] = integrate(C_temp,(z,0,1)) / bandwidth

plt.plot(r1[0,0,:], ugrid)
plt.show()
np.savetxt("filtered_solution_scheme=euler_z0=0.50_N=100000000_h=2.0e-05.txt", np.reshape(r1, (len(dt_list), len(t_obs_list)*len(ugrid))), delimiter="\t")
