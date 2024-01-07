using Plots, DelimitedFiles, Printf
using CUDA, Adapt, Test, BenchmarkTools, SpecialFunctions, LinearAlgebra

function new_Kernel_Estimator(x, Locparticle, bandwidth)
    # kernel methods
    function Gauss_kernel(u; d=1)
        return (2*pi)^(-0.5*d) * exp(-0.5*u*u)
    end

    function Epa_kernel(u; d=1)
        nu_d = 2*pi^(0.5*d) / d / gamma(0.5*d) 
        if u*u >= 1
            return 0.
        else
            return 0.5 / nu_d * (d+2) * (1-u*u)
        end
    end

    #Kernel = Gauss_kernel
    Kernel = Epa_kernel
    #cgrid = CUDA.zeros(length(xgrid))
    d = 1
    p = 0 

    if x<=bandwidth
        #print("reflection applies.")
        u = (Locparticle - x)/bandwidth
        p += Kernel(u)
        # the kernel span outside the boundary; max when the sample is at the boundary
        # aka the deflection length
        d_reflection = bandwidth - (Locparticle-0) 
        if x<d_reflection
            u_reflection = (Locparticle + x)/bandwidth
            p += Kernel(u_reflection)
        end
    elseif 1-x <= bandwidth
        #print("reflection applies.")
        u = (Locparticle - x)/bandwidth
        p += Kernel(u)
        d_reflection = bandwidth-(1 -Locparticle)
        if 1-x<d_reflection
            u_reflection = (2-x - Locparticle)/bandwidth
            p += Kernel(u_reflection)
        end
    else
        u = (Locparticle - x)/bandwidth
        p += Kernel(u)
    end
    return p
end

function new_gpu_compute_concentration_by_kernel(t,x_grid,Locparticles)
    Nparticles = length(Locparticles)
    x_avg = sum(Locparticles)/length(Locparticles)
    x_delta = Locparticles .- x_avg
    std_sample = sqrt(1/(Nparticles - 1)* dot(x_delta, x_delta))
    bandwidth = std_sample * Nparticles^(-0.2)
    @printf("Calculation of std_sample completed. x_avg=%.4f, minx=%.8f, maxx=%.8f, std_sample=%.4f.\n",x_avg, minimum(Locparticles), maximum(Locparticles), std_sample)
    @printf("The bandwdith is %.4f\n",bandwidth)

    # use the kernel method to compute the PDF
    pdf_grid = similar(x_grid)
    # copy data from cpu to gpu
    cu_bandwidth = adapt(CuArray, [bandwidth])
    cu_x_grid = adapt(CuArray, x_grid)
    cu_pdf_grid = adapt(CuArray, pdf_grid)
    cu_Locparticles = adapt(CuArray, Locparticles)
    numthreads = 1024  # This can be changed to length of a 3D grid if necessary.
    numblocks = ceil(Int, Nparticles/numthreads)               # temporarily

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
                # there are failed attemps dealing with atomic operations. 
                # Now it is suggested to use a low-level interface with pointer, and convert the numerical value to supporting types
                #CUDA.@cushow(Kernel_Estimator(x_grid[i], Locparticles[j], bandwidth[1]))
                #@inbounds CUDA.@atomic pdf_grid[i] += Kernel_Estimator(x_grid[i], Locparticles[j], bandwidth[1]) # this does not require the use of pointers
                
                # for some unknown reasons, the kernel returns NaN from time to time at a really low frequency
                # to avoid its inference, rule out these cases
                temp = Float64(Kernel_Estimator(x_grid[i], Locparticles[j], bandwidth[1]))
                #if !isnan(temp)
                CUDA.atomic_add!(pointer(pdf_grid,Int(i)), Float64(temp))
                #end
                #@inbounds CUDA.atomic_add!(pointer(pdf_grid, Int(i)), Float64(Kernel_Estimator(x_grid[i], Locparticles[j], bandwidth[1])))
                #@inbounds CUDA.@atomic pdf_grid[i] +=  Float64(Kernel_Estimator(x_grid[i], Locparticles[j], bandwidth[1]))
            end
        end
        return nothing
    end

    CUDA.@sync begin
       @cuda threads=numthreads blocks=numblocks gpu_kernel_estimator!(new_Kernel_Estimator,cu_pdf_grid, cu_x_grid, cu_Locparticles, cu_bandwidth)
    end
    return bandwidth, Array(cu_pdf_grid) ./ (Nparticles*bandwidth^1) # dimension d=1
end


dirname = "case_gpu_testing//diffusion//temp//outputbw_231225"
z0=0.5
dt=1e-4
N=10000
h=2.0e-5
scheme="RK4"
fname=@sprintf("%s//pLocation_output_scheme=%s_z0=%.2f_N=%d_h=%.1e_dt=%.1e.txt", dirname, scheme, z0, N, h,dt)
pLocation = readdlm(fname,'\t',Float64,'\n')
data = pLocation[:,1]
data_len = length(data)
nan_location= []
xgrid = 0:0.02:1.0
grid_len = length(xgrid)
pdf_cpu = zeros(grid_len)
bw = 0.0389763092730652


for i=1:data_len
    for j=1:grid_len
        pdf_cpu[j] += new_Kernel_Estimator(xgrid[j], data[i], bw)
    end
end

pdf_cpu = pdf_cpu ./ (data_len*bw)

temp = zeros(data_len)
for i=1:data_len
    temp[i] = new_Kernel_Estimator(xgrid[15], data[i], bw)
end

λ, pdf = new_gpu_compute_concentration_by_kernel(0,xgrid, data)

for i=1:grid_len
    if isnan(pdf[i])
        push!(nan_location, [i, xgrid[i]])
    end
end

display(nan_location)
