using CUDA

function kernel(a,b)
    index = blockDim().x * (blockIdx().x - 1) + threadIdx().x
    stride = blockDim().x * gridDim().x
    for j = index:stride:length(a)
        temp =  j % 10
        CUDA.atomic_add!(pointer(a,Int(temp)), Float64(1))
    end
    #CUDA.@atomic a[id] += id
    return nothing
end

l = 100000
a = CUDA.zeros(Float64, l)
threadnum = 1024
blocknum = cld(l, threadnum)
b = CuArray([100,100])
@cuda threads=1024 blocks=blocknum kernel(a,b)

println(a[1:20])
println(Array(b))