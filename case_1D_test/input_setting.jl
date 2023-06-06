using Pkg
Pkg.instantiate()
Pkg.activate("..")
using Particles
using Random
randpool = MersenneTwister(0) #Generate same random numbers every time

# define the geometry of the computation domain
H = 1

# initialize the dict d that stores all the computation setting
d = default_userdata()
print(d)

# set appropriate times
d["dt"] = 1e-3     #time-step
d["tstart"] = 0.0 
d["tend"] = 1.0

# set particles
N = 1000
d["nparticles"] = N

# all variables for one particle are collected in a vector
variables = ["z", "age"]
d["variables"] = variables

# define the background velocity field
function w(z)
    0.0
end
  