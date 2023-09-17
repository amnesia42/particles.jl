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
dt = 3e-5; z0=0.5; N = 1000000
fname=@sprintf("./diffusion/z0=%.2f_N=%d_dt=%.2e.txt", z0, N, dt)
#fname="end_distribution_pure_diffusion_dt=3.0e-05_tend=0.200000_Nparticles=100000.txt"
# fname = blablabla #manually copy and paste
data = readdlm(fname, '\t', Float64, '\n')
layout = @layout([a  b])

p = plot()
t_obs = [0.036, 0.072, 0.108, 0.1728]
pdf_times = zeros(length(xs), length(t_obs))
Nparticles = size(data, 1)
if kernel_type=="box"
    dx = 0.02
    xs = 0:dx:1.0
    for i = 1:length(t_obs)
        ti = t_obs[i]
        d = data[:, i]
        h, bin_edges = np.histogram(d, bins=xs, density=true)
        plot!(p, h, bin_edges[1:end-1] .+ dx/2, linewidth=2, thickness_scaling=1,label="t=$ti", size=(600,600),dpi=300)
    end
else
    dx = 0.02
    xs = 0:dx:1.0
    for i = 1:length(t_obs)
        ti = t_obs[i]
        d = data[:, i]
        #pdf_times[:, i] = compute_concentration_by_kernel(ti, xs, d; Kernel_type=kernel_type)
        pdf_times[:, i] = gpu_compute_concentration_by_kernel(ti,xs,d)
        plot!(pdf_times[:, i], xs, label="t=$(ti)s", linewidth=2, thickness_scaling=1, size=(600,600), dpi=300)
    end

end
xlabel!("concentration")
ylabel!("location x")
title!("kernel=$(kernel_type)")
figname=@sprintf("./diffusion/kernel=%s_z0=%.2f_N=%d_dt=%.2e.png", kernel_type,z0, N, dt)
savefig(p, figname)