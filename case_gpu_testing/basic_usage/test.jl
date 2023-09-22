# This program checks whether a GPU and Julia Packages exist.
using CUDA
using Test

function gpu_add1!(a,b)
    index = threadIdx().x
    stride = blockDim().x
    for i = index:stride:length(b)
        @inbounds a[i] += b[i]
    end
    return nothing
end


# prepare the data for computation
N = 1024
a = CUDA.fill(1.0f0,N)
b = CUDA.fill(2.0f0,N)

# kernel programming
@cuda gpu_add1!(a,b)
@test all(Array(a) .== 3.0f0)

# array programming
t = sum(a .* b)