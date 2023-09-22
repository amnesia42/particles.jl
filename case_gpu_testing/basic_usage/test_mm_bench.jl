# This program involves several matrix-matrix product and kernel
# and benchmarks the computation time.
using CUDA
using Test
using BenchmarkTools
using LinearAlgebra

function get_L1err(A,B)
    return println("L1error=",sum(abs.(A .- B))/length(A))
end

function mat_mul!(A,B,C)
    m = size(A,1)
    l = size(A,2)
    n = size(B,2)
    for i=1:m
        for j=1:n
            C[i,j] = 0
            for k=1:l
                C[i,j] += A[i,k]*B[k,j]
            end
        end
    end
    return nothing
end

function mat_mul_vec!(A,B,C)
    m = size(A,1)
    l = size(A,2)
    n = size(B,2)
    for i=1:m
        for j=1:n
            C[i,j] = dot(A[i,:], B[:, j])
        end
    end
    return nothing
end

function gpu_mm!(A,B,C)
    # matrix a and b must be compataible for multiplication
    row1,col1 = CUDA.size(A)
    #if !(col1 != row2
    #    CUDA.throw_api_error("The dimension of two matrices don't match.")
    #    return nothing
    #end
    index = blockDim().x * (blockIdx().x - 1) + threadIdx().x
    stride = blockDim().x * gridDim().x 
    for k=index:stride:length(C)
        # change linear indices toCartesian indices
        j = div(k, row1, RoundUp)
        i = k - (j-1)*row1
        for l=1:col1 
            @inbounds C[k] += A[i, l] * B[l, j]
        end
    end
    return nothing
end

function bench_mm!(A,B,C, numthreads)
    numblocks=ceil(Int, length(A)/numthreads)
    CUDA.@sync begin
        @cuda threads=numthreads blocks=numblocks gpu_mm!(A,B,C)
    end
end

function bench2_mm!(A,B,C)
    kernel = @cuda launch=false gpu_mm!(A,B,C)
    config = launch_configuration(kernel.fun)
    threads = min(N, config.threads)
    blocks = cld(N, threads)
    #println("blockDim=$(threads), numblocks=$(blocks)")
    CUDA.@sync begin
        kernel(A,B,C; threads, blocks)
    end
end

function bench3_mm!(A,B,C)
    CUDA.@sync begin
        @cuda gpu_mm!(A,B,C)
    end
end

# prepare the data for computation
N = 1000; row1=N; row2=N; col1=N; col2=N;
A = rand(Float64, N,N)
B = rand(Float64, N,N)
C_seq = fill(0., (N,N))
# preallocate an array for storing the result
A_d = CuArray(A)
B_d = CuArray(B)
C_d = CUDA.fill(0., (N,N)) 
C_array = CUDA.fill(0., (N,N))


# computation 
println("===========Benchmarking starts================")

# analytical result
println("Sequential library:")
@btime C = A*B

# array programming
println("Parallel - array programming:")
@btime CUDA.@sync C_array = A_d * B_d 

# kernel programming
# This gives similar result as above, but the test might be susceptible to floating-point error.
println("User defined gpu grid:")
numthreads = 512
println("N=$(N), blockDim=$(numthreads), numblocks=$(ceil(Int, N^2/numthreads))")
@btime bench_mm!(A_d, B_d, C_d, numthreads)


println("Program-defined gpu grid:")
C_d .= 0
kernel = @cuda launch=false gpu_mm!(A_d,B_d,C_d)
config = launch_configuration(kernel.fun)
threads = min(N, config.threads)
blocks = cld(N, threads)
println("N=$(N), blockDim=$(threads), numblocks=$(blocks)")
@btime bench2_mm!(A_d,B_d,C_d)


println("===========Benchmarking ends================")
