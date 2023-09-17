include("ConcentrationCalculationLibrary.jl")
using Plots, DelimitedFiles, Printf

function gpukernel_time_stepping(p, dt; h=0.01)
    index = blockDim().x * (blockIdx().x - 1) + threadIdx().x
    stride = blockDim().x * gridDim().x     
    for k=index:stride:length(p)
        @inbounds p[k] = get_next_location(p[index], dt[1]; h=h)
    end
    return nothing
end

# init
N = 1000000
z0 = 0.5
Tend = 0.2
#dt = 1e-3 
dt_list = [3e-4,1e-4,3e-5,1e-5,3e-6]
scheme="euler"
for dt in dt_list
    t_obs = [0.036, 0.072, 0.108, 0.1728]
    pLocation_output = z0 * ones(N,length(t_obs))

    # copy data from cpu to gpu
    numthreads = 1024
    numblocks = cld(N, numthreads)
    cu_pLocation = CuArray(z0*ones(N))
    cu_dt = CUDA.fill(dt)


    # serial code for computing trajectory
    CUDA.@time begin
        closestindex = 1
        for t=0:dt:Tend
            CUDA.@sync begin
                @cuda threads=numthreads blocks=numblocks gpukernel_time_stepping(cu_pLocation, cu_dt)
            end
            if closestindex<=4 && abs(t-t_obs[closestindex])<1e-10
                pLocation_output[:,closestindex] = Array(cu_pLocation) 
                global closestindex += 1
            end
        end
    end


    # output to result
    fname=@sprintf("./diffusion/error/particles_scheme=%s_z0=%.2f_N=%d_dt=%.2e.txt", scheme, z0, N, dt)
    open(fname, "a") do io
        writedlm(io, pLocation_output)
    end
end