using Plots, DelimitedFiles, Printf
t_obs = [0.036, 0.072, 0.108, 0.144]
dt_list = [3e-3,1e-3,3e-4, 1e-4, 3e-5]
fname=@sprintf("./diffusion/error/particles_scheme=%s_z0=%.2f_N=%d_Nrep=%d.txt", scheme, z0, N, Nrep)
d = readdlm(fname, '\t', Float64, '\n')


