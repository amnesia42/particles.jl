using CUDA

function kernel(a)
    id = threadIdx().x % 10 + 1
    CUDA.atomic_add!(pointer(a,id), 1.0)
    return nothing
end

l = 10
a = CUDA.zeros(Float64, l)

@cuda threads=100 kernel(a)

println(sum(Array(a)))