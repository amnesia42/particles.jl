# There are two ways of parallelizing. 
# The first way is to assign the computation of each grid point to different sm (bottleneck-grid points)
# The second way is to assign different data to different procesors. (bottleneck-data)
# Both work for a simple 1D example.
# The bottleneck for this problem will be the grid points, because 100k particles is already a lot.
# Try to wrap this into a function that inputs a particle distribution and outs a concentration approximation.
kernel_type = "Epa" # "box", "Epa", "Gauss" if changed to "Gauss", remember to change in the Kernel_Estimator function in ConcentrationCalculationLibrary.jl
                    # This is a temporary solution because I don't know how to pass a static reference to the device memory. 

include("ConcentrationCalculationLibrary.jl")

function gpu_compute_concentration_by_kernel(t,x_grid,Locparticles)
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
    return Array(cu_pdf_grid)
end

# concentration estimation parameter specification
z0=0.5; N = 1000000; scheme="euler"

t_obs = [0.036, 0.072, 0.108, 0.1728]
dt_list = [3e-4,1e-4,3e-5,1e-5,3e-6]
error_list = zeros(length(dt_list), length(t_obs))
dx = 0.02
xs = 0:dx:1.0
pdf_times = zeros(length(xs), length(t_obs))
analpdf_times = zeros(length(xs), length(t_obs))
Nparticles = size(data, 1)

# compute analytical solution
for j = 1:length(t_obs)
    t = t_obs[j]
    for i=1:length(xs)
        analpdf_times[i,j] = 0.
        for n=0:10
            δC= (2*n+1)*Pl(2*xs[i]-1, n)*Pl(2*z0-1,n)*exp(-6*n*(n+1)*t)
            analpdf_times[i,j] += δC
        end
    end
end

for k=1:length(dt_list)
    dt = dt_list[k]
    fname=@sprintf("./diffusion/error/particles_scheme=%s_z0=%.2f_N=%d_dt=%.2e.txt", scheme, z0, N, dt)
    data = readdlm(fname, '\t', Float64, '\n')

    for i = 1:length(t_obs)
        ti = t_obs[i]
        d = data[:, i]
        pdf_times[:, i] = gpu_compute_concentration_by_kernel(ti,xs,d)
    end

    error_list[k,:] = sum(abs, pdf_times-analpdf_times, dims=1) ./ length(xs)
end

# output to result
fname=@sprintf("./diffusion/error/L1error_scheme=%s_z0=%.2f_N=%d.txt", scheme, z0, N)
open(fname, "w") do io
    writedlm(io, error_list)
end

