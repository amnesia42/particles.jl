`julia> include("test_mm_bench.jl")
===========Benchmarking starts================
Sequential library:
  20.670 ms (2 allocations: 7.63 MiB)
Parallel - array programming:
  18.427 ms (91 allocations: 4.33 KiB)
User defined gpu grid:
N=1000, blockDim=512, numblocks=1954
  45.031 ms (55 allocations: 3.66 KiB)
Program-defined gpu grid:
N=1000, blockDim=768, numblocks=2
  240.182 ms (78 allocations: 4.69 KiB)
===========Benchmarking ends===============
`

`julia> include("test_mm.jl")
===========Computation starts================
Sequential library:
  0.023889 seconds (2 allocations: 7.629 MiB)
Parallel - array programming:
  0.145744 seconds (91 allocations: 4.328 KiB)
L1error=8.325039857481898e-14
false
User defined gpu grid:
N=1000, blockDim=512, numblocks=1954
  0.245857 seconds (56.21 k allocations: 4.670 MiB, 7.48% compilation time)
L1error=1.5860555890867544e-13
Program-defined gpu grid:
N=1000, blockDim=768, numblocks=2
  0.289376 seconds (4.90 k allocations: 350.833 KiB, 4.73% compilation time)
L1error=8.325039857481898e-14
User implemented sequnetial mm-product:
  0.528058 seconds (8.31 k allocations: 574.927 KiB, 3.25% compilation time)
false
L1error=1.586171833878325e-13
User implemented partly vectorized mm-product:
  2.295886 seconds (2.01 M allocations: 15.140 GiB, 20.21% gc time, 1.35% compilation time)
L1error=2.948894461951568e-14
===========Computation ends================`