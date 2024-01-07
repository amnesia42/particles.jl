import numpy as np

p_tilde = np.loadtxt("pdfoutput_kernel=Epa_scheme=euler_z0=0.50_N=100000000.txt", delimiter="\t")
print(p_tilde.shape)