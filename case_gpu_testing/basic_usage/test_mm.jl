# This program involves several matrix-matrix product and kernel
# and checks its accuracy.

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
println("===========Computation starts================")

# analytical result
println("Sequential library:")
@time C = A*B

# array programming
println("Parallel - array programming:")
@time CUDA.@sync C_array = A_d * B_d 
get_L1err(Array(C_array), C)
println(all(broadcast(==, Array(C_array), C))) #@test all(Array(C_d) .== C)

# kernel programming
# This gives similar result as above, but the test might be susceptible to floating-point error.
println("User defined gpu grid:")
numthreads = 512
println("N=$(N), blockDim=$(numthreads), numblocks=$(ceil(Int, N^2/numthreads))")
@time bench_mm!(A_d, B_d, C_d, numthreads)
get_L1err(Array(C_d), C)


println("Program-defined gpu grid:")
C_d .= 0
kernel = @cuda launch=false gpu_mm!(A_d,B_d,C_d)
config = launch_configuration(kernel.fun)
threads = min(N, config.threads)
blocks = cld(N, threads)
println("N=$(N), blockDim=$(threads), numblocks=$(blocks)")
@time bench2_mm!(A_d,B_d,C_d)
get_L1err(Array(C_array), C)

# Sequential mm-product implementation
println("User implemented sequnetial mm-product:")
@time mat_mul!(A,B,C_seq)
println(all(broadcast(==, C_seq, C))) 
get_L1err(C, C_seq)
println("User implemented partly vectorized mm-product:")
C_seq = fill(0., (N,N))
@time mat_mul_vec!(A,B,C_seq)
get_L1err(C, C_seq)
println("===========Computation ends================")

# println("===========Error test 1 starts================")
# N = 300
# println("N=$N")
# println("The following computation is slow because it doesn't explicitly specify numthreads and numblocks.")

# @time CUDA.@sync @cuda gpu_mm!(A_d,B_d,C_d)
# kernel = @cuda launch=false gpu_mm!(A_d,B_d,C_d)
# config = launch_configuration(kernel.fun)
# threads = min(N, config.threads)
# blocks = cld(N, threads)
# println("N=$(N), blockDim=$(threads), numblocks=$(blocks)")
# @time bench3_mm!(A_d,B_d,C_d)

# println("===========Error test 2 ends================")




# Deprecated codes


#CUDA.@time C_array = A_d * B_d 
#println(all(broadcast(==, Array(C_array), C))) #@test all(Array(C_d) .== C)
#C_result = Array(C_array)
#L1error = sum(abs.(broadcast(-, C_result, C)))
#println(L1error/N^2)


# kernel programming
# This gives similar result as above, but the test might be susceptible to floating-point error.
# @cuda threads=256 gpu_mm!(A_d, B_d, C_d)
# So we compute the L1 error instead.
# C_result = Array(C_d)
# println(size(C))
# println(size(C_result))
# L1error = sum(abs.(broadcast(-, C_result, C)))
# println(L1error/N)



# println("User defined gpu grid:")
# numthreads = 512
# println("N=$(N), blockDim=$(numthreads), numblocks=$(ceil(Int, N^2/numthreads))")
# @btime bench_mm!(A_d, B_d, C_d, numthreads)

# println("Program-defined gpu grid:")
# C_d .= 0
# kernel = @cuda launch=false gpu_mm!(A_d,B_d,C_d)
# config = launch_configuration(kernel.fun)
# threads = min(N, config.threads)
# blocks = cld(N, threads)
# println("N=$(N), blockDim=$(threads), numblocks=$(blocks)")
# @btime bench2_mm!(A_d,B_d,C_d)

# Benchmarking


# # Sequential mm-product implementation
# println("User implemented sequnetial mm-product:")
# @time mat_mul!(A,B,C_seq)
# println(C_seq)
# println(all(broadcast(==, C_seq, C))) 
# L1error = sum(abs.(C_seq .- C))
# println(L1error/N^2)
# println("User implemented partly vectorized mm-product:")
# C_seq = fill(0., (N,N))
# @time mat_mul_vec!(A,B,C_seq)
# println(C_seq)
# println(get_L1err(C, C_seq))