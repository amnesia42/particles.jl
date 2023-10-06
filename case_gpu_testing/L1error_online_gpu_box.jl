using Plots, DelimitedFiles, Printf
using LegendrePolynomials
# There are two ways of parallelizing. 
# The first way is to assign the computation of each grid point to different sm (bottleneck-grid points)
# The second way is to assign different data to different procesors. (bottleneck-data)
# Both work for a simple 1D example.
# The bottleneck for this problem will be the grid points, because 100k particles is already a lot.
# Try to wrap this into a function that inputs a particle distribution and outs a concentration approximation.
kernel_type = "box" # "box", "Epa", "Gauss" if changed to "Gauss", remember to change in the Kernel_Estimator function in ConcentrationCalculationLibrary.jl
                    # This is a temporary solution because I don't know how to pass a static reference to the device memory. 
scheme = "euler" # "euler", "m1", "heun","RK4" remember to change in the get_next_location function in ConcentrationCalculationLibrary.jl
# This is a temporary solution because I don't know how to pass a static reference to the device memory. 

include("ConcentrationCalculationLibrary.jl")
function gpukernel_time_stepping(p, dt; h=0.002)
    index = blockDim().x * (blockIdx().x - 1) + threadIdx().x
    stride = blockDim().x * gridDim().x     
    for k=index:stride:length(p)
        @inbounds p[k] = get_next_location(p[index], dt[1]; h=h)
    end
    return nothing
end

function get_particles_snapshot(z0, N, scheme, dt, Tend, t_obs)
    N_tobs = length(t_obs)
    pLocation_output = z0 * ones(N,N_tobs)

    # copy data from cpu to gpu
    numthreads = 1024
    numblocks = cld(N, numthreads)
    cu_pLocation = CuArray(z0*ones(N))
    cu_dt = CUDA.fill(dt)

    # serial code for computing trajectory
    closestindex = 1
    for t=0:dt:Tend
        CUDA.@sync begin
            @cuda threads=numthreads blocks=numblocks gpukernel_time_stepping(cu_pLocation, cu_dt)
        end
        if closestindex<=N_tobs && abs(t-t_obs[closestindex])<1e-10
            pLocation_output[:,closestindex] = Array(cu_pLocation) 
            closestindex += 1
        end
    end
    return pLocation_output
end

# concentration estimation parameter specification
z0=0.5; N = 100000000;Tend = 0.216;scheme="euler"
t_obs = [0.036, 0.072, 0.108, 0.144, 0.180, 0.216]
#t_obs = 3e-3:3e-3:0.216
dt_list = [3e-3,1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]
#t_obs = [0.036, 0.072, 0.108, 0.1728]

#dt_list = [3e-4,1e-4,3e-5,1e-5,3e-6]
dx = 0.02
xs = 0:dx:1.0
xs_true = (0+dx/2.0):dx:(1.0-dx/2.0)
error_output = zeros(length(dt_list), length(t_obs))
bandwidth_output = zeros(length(dt_list), length(t_obs))
pdf_output = zeros(length(dt_list), length(xs_true),length(t_obs))
analpdf_times = zeros(length(xs_true), length(t_obs))
Nparticles = N

# compute analytical solution
for j = 1:length(t_obs)
    t = t_obs[j]
    for i=1:length(xs_true)
        analpdf_times[i,j] = 0.
        for n=0:50
            δC= (2*n+1)*Pl(2*xs_true[i]-1, n)*Pl(2*z0-1,n)*exp(-6*n*(n+1)*t)
            analpdf_times[i,j] += δC
        end
    end
end

for i=1:length(dt_list)
    dt=dt_list[i]
    pdf_times = zeros(length(xs_true), length(t_obs))
    data = get_particles_snapshot(z0, N, scheme, dt, Tend, t_obs)
    for k = 1:length(t_obs)
        ti = t_obs[k]
        d = data[:, k]
        pdf_times[:, k] += Box_Estimator(xs,d)
        bandwidth_output[i,k] = get_bandwidth(d)
    end
    println("dt=",dt," computation finiished.")
    pdf_output[i,:,:] = pdf_times
    error_output[i,:] = sum(abs, pdf_times - analpdf_times, dims=1)
end

# for i=1:length(dt_list)
#     for j=1:length(t_obs)
#         for k=1:length(xs)
#             if k==1 || k==length(xs)
#                 error_output[j,i] += 0.5*error_list[k,j,i]
#             else
#                 error_output[j,i] += error_list[k,j,i]
#             end
#         end
#     end
# end


# output to result
fname=@sprintf("./diffusion/error/L1error_kernel=%s_scheme=%s_z0=%.2f_N=%d.txt", kernel_type,scheme, z0, N)
open(fname, "w") do io
    writedlm(io, error_output)
end

# output to result
fname=@sprintf("./diffusion/error/pdfoutput_kernel=%s_scheme=%s_z0=%.2f_N=%d.txt", kernel_type, scheme, z0, N)
open(fname, "w") do io
    writedlm(io, pdf_output)
end

# output to result
fname=@sprintf("./diffusion/error/bandwidthoutput_kernel=%s_scheme=%s_z0=%.2f_N=%d.txt", kernel_type, scheme, z0, N)
open(fname, "w") do io
    writedlm(io, bandwidth_output)
end
