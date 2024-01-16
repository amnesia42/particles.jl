import numpy as np
import matplotlib.pyplot as plt
#from math import Rational
from sympy import legendre, diff, integrate, factorial, exp, Rational,lambdify, Piecewise,simplify
from sympy.abc import x, t, u, z

def check_solution(C):
    C = simplify(C)
    C_sol = lambdify([u,t], C, "numpy") # make it act as a normal function
    print(C_sol(0.5,0.036))
    N_zspacing = 50
    zgrid = np.linspace(0,1,N_zspacing+1)
    #fig,ax = plt.subplots()
    for time in [0.036]:#,0.072#,0.108,0.144,0.180,0.216]:
        res = np.zeros((N_zspacing+1,))
        for i in range(N_zspacing+1):
            zloc = zgrid[i]
            res[i] =  C_sol(zloc, time) # C.evalf(subs={x:zloc, t:0.036})
        plt.plot(res, zgrid, label="t={:.3f}".format(time))
        print(res)
        plt.show()    
    return None

def check_solution_old(C):
    N_zspacing = 50
    zgrid = np.linspace(0,1,N_zspacing+1)
    res = np.zeros((N_zspacing+1,))
    for time in [0.036,0.072,0.108,0.144,0.180,0.216]:
        for i in range(N_zspacing+1):
            zloc = zgrid[i]
            res[i] =  C.evalf(subs={u:zloc, t:time})
        plt.plot(res, zgrid, label="t={:.3f}".format(time))
    plt.legend()
    plt.xlabel("Concentration C")
    plt.ylabel("Location z")
    plt.show()
    return None

################ SIMPLE TEST STARTS ###########################
print("SIMPLE TEST STARTS")
expr = 1 
print(type(expr))
for i in range(1,11):
    expr += x**i/factorial(i)
print(expr)
print(expr.evalf())
print(expr.evalf(subs={x: 1.0}))
print("SIMPLE TEST ENDS")
################ SIMPLE TEST ENDS ###########################
################ SOLUTION PLOT STARTS ###########################
print("SOLUTION PLOT STARTS")
N_legendre = 50 # number of legendre polynomial terms
K = Piecewise( (Rational(3,4)*(1 - x**2), x*x<=1), (0, True) )
C = 1
for n in range(1,N_legendre): 
    C += (2*n+1) * legendre(n, 2*u-1) * legendre(n,0) * exp(-6*n*(n+1)*t)
check_solution_old(C)
print("SOLUTION PLOT ENDS")
################ SOLUTION PLOT ENDS ###########################
################ INTEGRATION STARTS ###########################
# print("INTEGRATION STARTS")
# ugrid = np.linspace(0,0.5,26)
# r1 = np.zeros_like(ugrid)
# bandwidth = 0.03
# K = K.subs(x,(u-z)/bandwidth)
# C_integrand = C * K
# for index in range(1):
#     grid = ugrid[index]
#     C_temp = C_integrand.subs({u: grid, t: 0.036})
#     r1[index] = integrate(C_temp,(z,0,1))
# plt.plot(r1, ugrid)
# plt.show()
# check_solution(C_int_sol)
# C_int_sol_derivative = diff(C_int_sol, u)
# C_int_sol_derivative_lambdify = lambdify([x,t], C_int_sol_derivative, "numpy")
# check_solution(C_int_sol_derivative)
# print("INTEGRATION ENDS")
################ INTEGRATION ENDS ###########################
