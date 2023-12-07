import numpy as np
import matplotlib.pyplot as plt
#from math import Rational
from sympy import legendre, diff, integrate, factorial, exp, Rational,lambdify
from sympy.abc import x, t, u

def check_solution(C):
    C_sol = lambdify([x,t], C, "numpy") # make it act as a normal function
    print(C_sol(0.5,0.036))
    N_zspacing = 50
    zgrid = np.linspace(0,1,N_zspacing+1)
    res = np.zeros((N_zspacing+1,))
    for i in range(N_zspacing+1):
        zloc = zgrid[i]
        res[i] =  C_sol(zloc,0.036) # C.evalf(subs={x:zloc, t:0.036})
    plt.plot(res, zgrid)
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
K = Rational(3,4)*(1 - x**2)
C = 1
for n in range(1,N_legendre): 
    C += (2*n+1) * legendre(n, 2*x-1) * legendre(n,0) * exp(-6*n*(n+1)*t)
check_solution(C)
print("SOLUTION PLOT ENDS")
################ SOLUTION PLOT ENDS ###########################
################ INTEGRATION STARTS ###########################
print("INTEGRATION STARTS")
bandwidth = 0.03
K = K.subs(x,(x-u)/bandwidth)
C_integrand = C * K
C_int_sol = integrate(C_integrand,(x,0,1))
check_solution(C_int_sol)
C_int_sol_derivative = diff(C_int_sol, u)
C_int_sol_derivative_lambdify = lambdify([x,t], C_int_sol_derivative, "numpy")
check_solution(C_int_sol)
print("INTEGRATION ENDS")
################ INTEGRATION ENDS ###########################
