using CUDA

# number of cores/SM 
# check website: eg. RTX 3050 has at least 2304 CUDA cores

println(attribute(device(),CUDA.DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X)) # DIM_X/DIM_Y:1024 DIM_Z:64

println("Max number of threads per block:",attribute(device(),CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK))       #1024
println("Max number of threads per sm:",attribute(device(),CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR)) #1536
println("Max number of blocks per sm:",attribute(device(),CUDA.DEVICE_ATTRIBUTE_MAX_BLOCK_PER_MULTIPROCESSOR))    #16
