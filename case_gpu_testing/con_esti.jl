# Try to wrap this into a function that inputs a particle distribution and outs a concentration approximation.

using CUDA, Test, BenchmarkTools
using LinearAlgebra, SpecialFunctions
using DelimitedFiles
using Plots
using PyCall, Printf
np = pyimport("numpy")

# diffusion functions not in a hurry
    # kernel methods
function Gauss_kernel(u; d=1)
    return (2*pi)^(-0.5*d) * exp(-0.5*u*u)
end

function Epa_kernel(u; d=1)
    nu_d = 2*pi^(0.5*d) / d / gamma(0.5*d) 
    if u*u > 1
        return 0.
    else
        return 0.5 / nu_d * (d+2) * (1-u*u)
    end
end
# kernel estimator using the reflection method to have order-1 consistency
function Kernel_Estimator(x, Locparticles, bandwidth; Kernel_type="Epa")
    if Kernel_type=="Epa"
       Kernel = Epa_kernel
    elseif Kernel_type=="Gauss"
       Kernel = Gauss_kernel
    end
    d = 1
    Nparticles = length(Locparticles)
    p = 0
    if x<bandwidth
        #print("reflection applies.")
        for i in eachindex(Locparticles)
            u = (Locparticles[i] - x)/bandwidth
            p += Kernel(u)
            # the kernel span outside the boundary; max when the sample is at the boundary
            # aka the deflection length
            d_reflection = bandwidth-(Locparticles[i]-0) 
            if x<d_reflection
                u_reflection = (Locparticles[i] + x)/bandwidth
                p+= Kernel(u_reflection)
            end
        end
        temp = p/(Nparticles*bandwidth^d)
    elseif 1-x< bandwidth
        #print("reflection applies.")
        for i in eachindex(Locparticles)
            u = (Locparticles[i] - x)/bandwidth
            p += Kernel(u)
            d_reflection = bandwidth-(1 -Locparticles[i])
            if 1-x<d_reflection
                u_reflection = (2-x - Locparticles[i])/bandwidth
                p+= Kernel(u_reflection)
            end
        end
        temp = p/(Nparticles*bandwidth^d)
    else
        for i in eachindex(Locparticles)
            u = (Locparticles[i] - x)/bandwidth
            p+= Kernel(u)
        end
        temp = p/(Nparticles*bandwidth^d)
    end
    return temp
end

function compute_concentration_by_kernel(t,x_grid, Locparticles; Kernel_type="Epa")
    Nparticles = length(Locparticles)
    x_avg = sum(Locparticles)/length(Locparticles)
    x_delta = Locparticles .- x_avg
    std_sample = sqrt(1/(Nparticles - 1)* dot(x_delta, x_delta))
    bandwidth = std_sample * Nparticles^(-0.2)
    @printf("Calculation of std_sample completed. std_sample=%.4f.\n",std_sample)
    @printf("The bandwdith is %.4f\n",bandwidth)

    # use the kernel method to compute the PDF
    pdf_grid = similar(x_grid)
    for i in eachindex(x_grid)
        pdf_grid[i] = Kernel_Estimator(x_grid[i], Locparticles, bandwidth; Kernel_type=Kernel_type)
    end
    @printf("Calculation of concentration profile at t=%.4f completed.\n", t)
    return pdf_grid
end

# concentration estimation parameter specification
kernel_type = "Epa" # "box", "Epa", "Gauss"
Nparticles::Int = 1000

# load the data
fhead = "end_distribution_pure_diffusion"
x0 = 0.5 # particles releasing location
dt = 3e-5   
t_end = 0.2
N_step = round(Int, t_end/dt) # number of times to update the particle location
t_obs = [0.036, 0.072, 0.108, 0.1728]
#t_obs = [1.0,2.0]
kthstep_obs = broadcast(x->round(Int, x/dt), t_obs) # note that no -1 is needed because this computation returns the result after the kth update
location_obs = zeros(Nparticles, length(t_obs))
#fname=string(fhead,"_","dt=$dt","_","tend=$t_end","_","Nparticles=$N_particles",".txt")
fname="distribution_pure_diffusion_dt=1.0e-04_tend=0.200000_Nparticles=1000.txt"
dname="F:\\Master_Thesis\\particles model\\1D advection diffusion\\parabola_result"
fname = joinpath(dname, fname)
# fname = blablabla #manually copy and paste
data = readdlm(fname, ' ', Float64, '\n')
layout = @layout([a  b])

#d = d[1:Nparticles,1]
# for plot_index in eachindex(Nparticles_list)
#     Nparticles_plot = Nparticles_list[plot_index]
#     if Nparticles_plot != Nparticles
#         location_obs = zeros(Nparticles_plot, length(t_obs))
#     end
#     if kernel_type=="box"
#         for ti in range(len(t_obs))
#             zs = np.linspace(0,1,101)
#             p
#         end
#     end
# end
p = plot()
dx = 0.01
xs = 0:dx:1.0
pdf_times = zeros(length(xs), length(t_obs))
if kernel_type=="box"
    for i = 1:length(t_obs)
        ti = t_obs[i]
        d = data[:, i]
        h, bin_edges = np.histogram(d, bins=xs, density=true)
        plot!(p, h, bin_edges[1:end-1] .+ dx/2, label="t=$ti")
    end
else
    for i = 1:length(t_obs)
        ti = t_obs[i]
        d = data[:, i]
        pdf_times[:, i] = compute_concentration_by_kernel(ti, xs, d; Kernel_type=kernel_type)
        plot!(pdf_times[:, i], xs, label="t=$(ti)s", size=(600,600))
    end

end
xlabel!("concentration")
ylabel!("location x")
title!("kernel=$(kernel_type)")
#savefig(p, "./case_gpu_tesing/figures/kernel=$(kernel_type)_Nparticles=$(Nparticles).png")
savefig(p, "./figures/kernel=$(kernel_type)_Nparticles=$(Nparticles).svg")