import numpy as np
import matplotlib.pyplot as plt
import random

mu = 0.5
sigma = 1.0
Nlist = [2**i for i in range(3,10)]
errorlist = np.zeros((len(Nlist),))
Tend = 1
M = 15000000
for i in range(len(Nlist)):
    N = Nlist[i]
    dt = Tend/N
    S = np.ones((M,))
    for j in range(N):
        dW = np.random.randn(M) * np.sqrt(dt)
        S = S + mu*S*dt + sigma*S*dW
    errorlist[i] = np.abs(np.mean(S) - np.exp(mu*Tend))
    print("N={:d}, err={:e}".format(int(N), np.abs(np.mean(S) - np.exp(mu*Tend))))
    print("X_T={:e}, X(T)={:e}".format(np.mean(S), np.exp(mu*Tend)))

for i in range(len(Nlist) - 1):
    tg = abs( (np.log10(errorlist[i+1]) - np.log10(errorlist[i]) ) / ( np.log10(Nlist[i+1]) - np.log10(Nlist[i]) )  )
    tg_text = "{:.3f}".format(tg)
    plt.text(0.5*(Nlist[i]+Nlist[i+1]), 0.5*(errorlist[i]+errorlist[i+1]), tg_text)

plt.loglog(Nlist, errorlist, marker='o')
plt.xlabel("N")
plt.ylabel("|EX_T-EX(T)|")
plt.savefig("./figures/pyweakeuler.png", dpi=300)