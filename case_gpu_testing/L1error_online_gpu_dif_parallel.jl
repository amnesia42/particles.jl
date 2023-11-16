using Plots, DelimitedFiles, Printf
using LegendrePolynomials
# There are two ways of parallelizing. 
# The first way is to assign the computation of each grid point to different sm (bottleneck-grid points)
# The second way is to assign different data to different procesors. (bottleneck-data)
# Both work for a simple 1D example.
# The bottleneck for this problem will be the grid points, because 100k particles is already a lot.
# Try to wrap this into a function that inputs a particle distribution and outs a concentration approximation.
kernel_type = "Epa" # "box", "Epa", "Gauss" if changed to "Gauss", remember to change in the Kernel_Estimator function in ConcentrationCalculationLibrary.jl
                    # This is a temporary solution because I don't know how to pass a static reference to the device memory. 
scheme = "euler" # "euler", "m1", "heun","RK4" remember to change in the get_next_location function in ConcentrationCalculationLibrary.jl
# This is a temporary solution because I don't know how to pass a static reference to the device memory. 

include("ConcentrationCalculationLibrary.jl")
function gpukernel_time_stepping(p, dt; h=0.0002)
    index = blockDim().x * (blockIdx().x - 1) + threadIdx().x
    stride = blockDim().x * gridDim().x     
    for k=index:stride:length(p)
        @inbounds p[k] = get_next_location(p[index], dt[1]; h=h)
    end
    return nothing
end

function gpu_compute_concentration_by_kernel(t,x_grid,Locparticles)
    Nparticles = length(Locparticles)
    x_avg = sum(Locparticles)/length(Locparticles)
    x_delta = Locparticles .- x_avg
    std_sample = sqrt(1/(Nparticles - 1)* dot(x_delta, x_delta))
    bandwidth = std_sample * Nparticles^(-0.2)
    @printf("Calculation of std_sample completed. std_sample=%.4f.\n",std_sample)
    @printf("The bandwidth is %.4f\n",bandwidth)

    # use the kernel method to compute the PDF
    pdf_grid = similar(x_grid)
    # copy data from cpu to gpu
    cu_bandwidth = CuArray([bandwidth])
    cu_x_grid = CuArray(x_grid)
    cu_pdf_grid = CuArray(pdf_grid)
    cu_Locparticles = CuArray(Locparticles)
    numthreads = length(x_grid)  # This can be changed to length of a 3D grid if necessary.
    numblocks = 1                # temporarily

    # check the Kernel Estimator can be passes a bitstype variable
    println(isa(Kernel_Estimator, Function))
    println(isbitstype(Kernel_Estimator))

    function gpu_kernel_estimator!(Kernel_Estimator, pdf_grid, x_grid,Locparticles, cu_bandwidth)
        index = blockDim().x * (blockIdx().x - 1) + threadIdx().x
        @inbounds pdf_grid[index] = Kernel_Estimator(x_grid[index], Locparticles, cu_bandwidth[1])
        return nothing
    end

    CUDA.@time CUDA.@sync begin
       @cuda threads=numthreads blocks=numblocks gpu_kernel_estimator!(Kernel_Estimator, cu_pdf_grid, cu_x_grid, cu_Locparticles, cu_bandwidth)
    end
    return bandwidth, Array(cu_pdf_grid)
end

function null_func()
    return nothing
end

function new_gpu_compute_concentration_by_kernel(t,x_grid,Locparticles)
    Nparticles = length(Locparticles)
    x_avg = sum(Locparticles)/length(Locparticles)
    x_delta = Locparticles .- x_avg
    std_sample = sqrt(1/(Nparticles - 1)* dot(x_delta, x_delta))
    bandwidth = std_sample * Nparticles^(-0.2)
    @printf("Calculation of std_sample completed. std_sample=%.4f.\n",std_sample)
    @printf("The bandwdith is %.4f\n",bandwidth)

    # use the kernel method to compute the PDF
    pdf_grid = similar(x_grid)
    # copy data from cpu to gpu
    cu_bandwidth = adapt(CuArray, [bandwidth])
    cu_x_grid = adapt(CuArray, x_grid)
    cu_pdf_grid = adapt(CuArray, pdf_grid)
    cu_Locparticles = adapt(CuArray, Locparticles)
    numthreads = 1024  # This can be changed to length of a 3D grid if necessary.
    numblocks = ceil(Int, Nparticles/1024)               # temporarily

    # check the Kernel Estimator can be passes a bitstype variable
    # println(isa(new_Kernel_Estimator, Function))
    # println(isbitstype(new_Kernel_Estimator))
    # println(isa(Kernel_Estimator, Function))
    # println(isbitstype(Kernel_Estimator))

    function gpu_kernel_estimator!(Kernel_Estimator, pdf_grid, x_grid, Locparticles, bandwidth)
        index = blockDim().x * (blockIdx().x - 1) + threadIdx().x
        stride = blockDim().x * gridDim().x
        for j = index:stride:length(Locparticles)
            for i = 1:length(pdf_grid)
                @inbounds CUDA.atomic_add!(pointer(pdf_grid,i), Kernel_Estimator(x_grid[i], Locparticles[j], bandwidth[1]))
            end
        end
        return nothing
    end

    CUDA.@sync begin
       @cuda threads=numthreads blocks=numblocks gpu_kernel_estimator!(new_Kernel_Estimator,cu_pdf_grid, cu_x_grid, cu_Locparticles, cu_bandwidth)
    end
    return bandwidth, Array(cu_pdf_grid) ./ (Nparticles*bandwidth^1) # dimension d=1
end

function get_particles_snapshot(z0, N, scheme, dt, Tend, t_obs)
    N_tobs = length(t_obs)
    pLocation_output = z0 * ones(N,N_tobs)

    # copy data from cpu to gpu
    numthreads = 1024
    numblocks = cld(N, numthreads)
    cu_pLocation = adapt(CuArray, z0*ones(N))
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
error_output = zeros(length(dt_list), length(t_obs))
bandwidth_output = zeros(length(dt_list), length(t_obs))
pdf_output = zeros(length(dt_list), length(xs),length(t_obs))
analpdf_times = zeros(length(xs), length(t_obs))
Nparticles = N

# compute analytical solution
for j = 1:length(t_obs)
    t = t_obs[j]
    for i=1:length(xs)
        analpdf_times[i,j] = 0.
        for n=0:50
            δC= (2*n+1)*Pl(2*xs[i]-1, n)*Pl(2*z0-1,n)*exp(-6*n*(n+1)*t)
            analpdf_times[i,j] += δC
        end
    end
end

for i=1:length(dt_list)
    dt=dt_list[i]
    pdf_times = zeros(length(xs), length(t_obs))
    @time data = get_particles_snapshot(z0, N, scheme, dt, Tend, t_obs)
    for k = 1:length(t_obs)
        ti = t_obs[k]
        d = data[:, k]
        bandwidth_output[i,k], pdf_times[:, k] = new_gpu_compute_concentration_by_kernel(ti,xs,d)
    end
    println("dt=",dt," computation finiished.")
    pdf_output[i,:,:] = pdf_times
    error_output[i,:] = sum(abs, pdf_times - analpdf_times, dims=1)
end


# output to result
fname=@sprintf("./diffusion/temp/L1error_scheme=%s_z0=%.2f_N=%d.txt", scheme, z0, N)
open(fname, "w") do io
    writedlm(io, error_output)
end

# output to result
fname=@sprintf("./diffusion/temp/pdfoutput_scheme=%s_z0=%.2f_N=%d.txt", scheme, z0, N)
open(fname, "w") do io
    writedlm(io, pdf_output)
end
