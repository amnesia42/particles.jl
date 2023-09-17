include("ConcentrationCalculationLibrary.jl")
using Plots, DelimitedFiles, Printf

# init
N = 5000
z0 = 0.5
Tend = 0.2
dt = 1e-4
t_obs = [0.036, 0.072, 0.108, 0.1728]
pLocation = z0 * ones(N)
pLocation_output = z0 * ones(N,length(t_obs))


# serial code for computing trajectory
closestindex = 1
for t=0:dt:Tend
    for i=1:N
        pLocation[i] = get_next_location(pLocation[i], dt)
    end
    if abs(t-t_obs[closestindex])<1e-10
        pLocation_output[:,closestindex] = pLocation
        if closestindex==4
            global closestindex = 1
        else 
            global closestindex += 1
        end
    end
end


# output to result
fname=@sprintf("./diffusion/z0=%.2f_N=%d_dt=%.2e.txt", z0, N, dt)
open(fname, "a") do io
    writedlm(io, pLocation_output)
end