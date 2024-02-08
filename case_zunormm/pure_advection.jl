# Initialize the notebook
using Pkg
#Pkg.activate("..")
#include("../src/Particles.jl")
using Particles
using Plots
using Zarr
using ZipFile
using DiskArrays
using Printf


# unzip the flow data (if not already done)

# function unzip(file,exdir="")
# extract all files from zip file to exdir
function unzip(file,exdir="")
    fileFullPath = isabspath(file) ?  file : joinpath(pwd(),file)
    basePath = dirname(fileFullPath)
    outPath = (exdir == "" ? basePath : (isabspath(exdir) ? exdir : joinpath(pwd(),exdir)))
    isdir(outPath) ? "" : mkdir(outPath)
    zarchive = ZipFile.Reader(fileFullPath)
    for f in zarchive.files
        fullFilePath = joinpath(outPath,f.name)
        if (endswith(f.name,"/") || endswith(f.name,"\\"))
            mkdir(fullFilePath)
        else
            write(fullFilePath, read(f))
        end
    end
    close(zarchive)
end

if !isdir("data/ZUNORMM_map_fullgrid_v3.zarr") 
    unzip("data/ZUNORMM_map_fullgrid_v3.zip","")
end

Zarr_data = ZarrData("data","ZUNORMM_map_fullgrid_v3.zarr")
t0=get_reftime(Zarr_data)
h=initialize_interpolation(Zarr_data,"waterlevel",t0)
depth = initialize_interpolation(Zarr_data,"waterdepth",t0)
u=initialize_interpolation(Zarr_data,"x_velocity",t0)
v=initialize_interpolation(Zarr_data,"y_velocity",t0)
w=initialize_interpolation(Zarr_data,"z_velocity_center", t0)
s=initialize_interpolation(Zarr_data,"salinity",t0)

# model parameters
x_left = 3.4 # in degrees
x_right = 5.6
y_left = 51.6
y_right = 52.5
z_surface = 0.0
z_bottom = -40.0
t_start=0.0
t_step=2*3600.0 # inseconds
#t_stop=15*24*3600.0 # end of dataset at 15 days
t_stop=12*3600.0 # can be 12 hours

# particle-tracking parameters
# Load default settings and adjust

# collected configuration is in Dict d
d = default_userdata() # start with some defaults
# all variables for one particle are collected in a vector
variables = ["x", "y", "z", "age"]
d["variables"] = variables
# initial position of the particles
m = length(variables)
# simulation time
d["dt"] = t_step/12     #time-step
d["tstart"] = t_start #start after 2 cycles
d["tend"]   = t_stop/3 #3600.0

# keep some output in memory
d["keep_particles"] = true #keep results in memory (bad idea for a large run)
d["keep_particle_times"] = collect(tstart:30.0:tend)

# personally set parameters
# the following options are to be added in the default_userdata() function
# if to apply the algorithm to prevent the algorithm from crossing the boundary or not
d["is_cross_bc"] = false # if added, the algorithm will terminate once the particle is out of the boundary
d["is_apply_cross_bc_test"] = false # if applied, only one particle is used and all virtual trajectories are also recorded
if d["is_apply_cross_bc_test"]
    d["the_virtual_trajectory"] = d["particles"]
end
# to be added according to the actual problem
d["is_particle_in_domain"] = x -> nothing
d["recursion_depth"] = 10
d["nparticles_bound"] = 200000
d["particles"] = zeros(length(variables), d["nparticles_bound"])
d["first_active_index"] = 1
d["first_inactive_index"] = 1
d["isactive_particles"] = falses(d["nparticles_bound"],)

# Randompool setup and restart
randpool = MersenneTwister(123) #Generate same random numbers every time
d["rng"] = randpool

# add the displacement update algorithm
include("displacement_update.jl")

# use the original grid in sampling and particle-tracking
# add the subroutine that checks whether the particle is inside the domain
# compulsory input             
function is_particle_in_domain(p)
    x,y,z,t = p
    is_inside_xdir = x_left ≤ x ≤ x_right
    is_inside_ydir = y_left ≤ y ≤ y_right
    is_inside_zdir = z_bottom ≤ z ≤ h(x,y_middle,z,t)
    is_apply_reflection = !is_inside_zdir            # in the tidal flume case, only particles going from above or below shall be reflected
    return is_inside_xdir && is_inside_zdir, is_apply_reflection 
end
d["is_particle_in_domain"] = is_particle_in_domain