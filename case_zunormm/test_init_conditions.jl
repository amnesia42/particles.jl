# Initialize the notebook
using Pkg
Pkg.activate("..")
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

if !isdir("./ZUNORMM_map_fullgrid_v2.zarr") 
    unzip("./ZUNORMM_map_fullgrid_v2.zip","")
end

# test case initialization 
# Load default settings and adjust
# collected configuration is in Dict d
d = default_userdata() # start with some defaults


Zarr_data = ZarrData(".","ZUNORMM_map_fullgrid_v2.zarr")
t0=get_reftime(Zarr_data)
h=initialize_interpolation(Zarr_data,"waterlevel",t0)
depth = initialize_interpolation(Zarr_data,"waterdepth",t0)
u=initialize_interpolation(Zarr_data,"x_velocity",t0)
v=initialize_interpolation(Zarr_data,"y_velocity",t0)
function w(x,y,z,t)
    return 0.0 # careful that x and y are recorded in 
end
s=initialize_interpolation(Zarr_data,"salinity",t0)


# find the x&y limit of the data
xy_bbox=Zarr_data.xy_grid.bbox
println("bbox=$(xy_bbox)") #bbox=[3.4773802907221345, 5.566420993890162, 51.66065856527758, 52.466444933890344]

# find t grid
t=Zarr_data.file.arrays["time"][:] #"seconds since 2022-04-01 00:00:00 +00:00"
d["tend"] = t[end]
d["dt"] = (t[2]-t[1])/ 6.0
println(t0) #"seconds since 2022-04-16 00:00:00 +00:00"
println( (t[2]-t[1])/3600.0 ) #2022-04-16 step=2hrs until 2022-05-01
println( (t[end]-t[2])/3600.0/24.0) #2022-04-16 duraton~=15hrs until 2022-05-01

# find x and y grid
# both are uniform, so store the grid size
d["x_center"]=Zarr_data.file.arrays["x_center"][:]
d["y_center"]=Zarr_data.file.arrays["y_center"][:]
d["dx"] = d["x_center"][2] - d["x_center"][1]
d["dy"] = d["y_center"][2] - d["y_center"][1]
d["z_iface_3d"] = Zarr_data.file.arrays["z_iface_3d"][:,:,:,:]

# # plot multiple z_center used in FVM computation
z_center3d = 0.01 .* Array(Zarr_data.file.arrays["z_center_3d"][:,:,:,1]) #999 x 499 x 46 x 181
z_center3d[z_center3d .>= 99.0] .= NaN

# # plot multiple z_surface used in FVM computation
z_iface3d = 0.01 .* Array(Zarr_data.file.arrays["z_iface_3d"][:,:,:,1]) #999 x 499 x 47 x 181
z_iface3d[z_iface3d .>= 99.0] .= NaN

# settings for this experiment
n = 0 # number of particles
d["nparticles"] = n

# all variables for one particle are collected in a vector
variables = ["x", "y", "z", "age"]
d["variables"] = variables
# initial position of the particles
m = length(variables)
p = zeros(m, n)
d["particles"] = p # initial values
# simulation time
#d["dt"] = 2*3600/6     # 2h/6, same as above
d["tstart"] = 0.0  #
d["tend"]   = 12*3600 # half a day
tstart=d["tstart"] 
tend=d["tend"]
# write to netcdf
d["write_maps_times"] = collect(0.0:d["dt"]*6:tend)
d["write_maps"] = false #do not write to netcdf
# write plots to file
d["plot_maps_times"] = collect(0.0:d["dt"]*6:tend)
d["plot_maps"] = false # do not make png figures
d["plot_maps_size"]=(900,300)

# keep some output in memory
d["keep_particles"] = true #keep results in memory (bad idea for a large run)
d["keep_particle_times"] = collect(tstart:d["dt"]*6:tend)

# add the function that drives movement of particles
# Here is the equation that we want to solve for the particles

"""
   !f(ds,s,t,i,d)

Dynamic model, computes as ds the function f in the equation ds=f(s,t)dt+g(s,t)dw 
for s at current time t for particle i and possibly using data/functions from d of type userdata.
"""
function f!(∂s, s, t, i, d)
    # constant defintion
    R_Earth = 6371000
    deg2rad = pi/180
    rad2deg = 180/pi
    # change velocity in m/s to degrees
    x, y, z, age = s
    # dx/dt=u
    ∂s.x = rad2deg * u(x, y, z, t) / (R_Earth*cos(deg2rad*y)) 
    # dy/dt=v 
    ∂s.y = rad2deg * v(x, y, z, t) / R_Earth
    # dz/dt=0
    ∂s.z = w(x, y, z, t)
    # age=(t-t0)
    ∂s.t = 1.0
end
d["f"] = f!

#
# TODO: add function g! to simulate diffusion to simulation routine
#

"""
   !g(ds,s,t,i,d)

   Dynamic model, computes as ds the function g in the equation ds=f(s,t)dt+g(s,t)dw 
   for s at current time t for particle i and possibly using data/functions from d of type userdata.
"""
function g!(∂s, s, t, i, d)
    # constant defintion
    R_Earth = 6371000
    deg2rad = pi/180
    rad2deg = 180/pi
    x, y, z, age = s
    # dx/dt=u
    ∂s.x = 0.0
    # dy/dt=v
    ∂s.y = 0.0
    # dz/dt=0
    ∂s.z = 0.0
    # age=(t-t0)
    ∂s.t = 0.0
end
d["g"] = g!

# Load data to check if particles is inside 
# the domain with the interpolated function

# load the z grid
z_center3d = 0.01 .* Array(Zarr_data.file.arrays["z_center_3d"][:,:,:,1]) #999 x 499 x 46 x 181
z_center3d[z_center3d .<= -99.0] .= NaN

# # plot multiple z_surface used in FVM computation
z_iface3d = 0.01 .* Array(Zarr_data.file.arrays["z_iface_3d"][:,:,:,1]) #999 x 499 x 47 x 181
z_iface3d[z_iface3d .<= -99.0] .= NaN

function is_particle_in_domain(p)
    x,y,z,t = p
    is_inside_xdir = xy_bbox[1] ≤ x ≤ xy_bbox[2]
    is_inside_ydir = xy_bbox[3] ≤ y ≤ xy_bbox[4]
    this_depth = depth(x,y,z,t)
    this_level = level(x,y,z,t)
    is_inside_zdir = (-1)*(this_depth-this_level) ≤ z ≤ this_level
    is_apply_reflection = !is_inside_zdir            # in the tidal flume case, only particles going from above or below shall be reflected
    return is_inside_xdir && is_inside_ydir && is_inside_zdir, is_apply_reflection 
end
d["is_particle_in_domain"] = is_particle_in_domain

function get_cell_size(t, xleft, xright, yleft, yright, zleft, zright, xcellnum, ycellnum, zcellnum)
    function get_cellcenters(left, right, gridlength)
        grid = range(left, right, length=gridlength)
        cellcenter = (grid[1:end-1] + grid[2:end]) ./ 2
        return cellcenter
    end
    xc = get_cellcenters(xleft, xright, xcellnum)
    yc = get_cellcenters(yleft, yright, ycellnum)
    zc = get_cellcenters(zleft, zright, zcellnum)
    wet_cellnumber = 0
    for x in xc
        for y in yc
            for z in zc
                if ~((s(x,y,z,t) - 9999) < 1e-10)
                    # if salinity exists
                    wet_cellnumber += 1
                end
            end
        end
    end
    # the appropritate cell size should be R*latitude*cos(\theta) * R*longtitude * z 
    # first simplifed version is latitude*cos(\theta) * longtitude * z 
    # second simplifed version is cos(\theta) * z (assume dx and dy are fixed)
    return wet_cellnumber/((xcellnum-1)*(ycellnum-1)*(zcellnum-1)) * cosd(0.5*(yright+yleft))*(zright-zleft)
end

function count_particles(xc_points, yc_points, z_ifaces, particles, isactivelist)
    # the input xc/yc must be compatible with z_ifaces
    # return concentration field defined in specified cell centers
    NumPar_field = zeros(length(xc_points), length(yc_points), length(z_ifaces)-2)
    active_particles_indexes = findall(isactivelist)
    for i in active_particles_indexes
        isindomain,_ = is_particle_in_domain(particles[:, i])
        if isindomain
            x, y, z, t = particles[:,i]
            box_xindex = findfirst(var->(x+d["dx"]/2)>=var, xc_points) # terminate searching at the first hit
            box_yindex = findfirst(var->(y+d["dy"]/2)>=var, yc_points) 
            box_zindex = findfirst(var->z>=var, z_ifaces[box_xindex,box_yindex,:]) # special data structure
            if !isnothing(box_xindex) 
                NumPar_field[box_xindex, box_yindex, box_zindex] += 1 # index match the upperbound, so minus 1 here
            end
        end
    end
    return NumPar_field
end
d["count_particles"] = count_particles

function approximate_concentration(t,xc_points, yc_points, z_ifaces, particles, isactivelist)
    NumPar_field = count_particles(xc_points, yc_points, z_ifaces, particles, isactivelist)
    cellsize = similar(NumPar_field)
    for i in eachindex(xc_points)
        x = xc_points[i]
        for j in eachindex(yc_points)
            y = yc_points[j]
            for k=1:length(z_ifaces)-1
                cellsize[i,j,k] = get_cell_size(t,x-d["dx"]/2,x+d["dx"]/2,
                                  y-d["dy"]/2,y+d["dy"]/2,z_ifaces[i,j,k],z_ifaces[i,j,k+1],
                                  10,10,10)
            end
        end
    end
    return NumPar_field ./ cellsize
end
d["approximate_concentration"] = approximate_concentration

function release_particles(d,p,t)
    #
end
d["release_particles"] = release_particles


function initialize_particles(d, p)    
    # check expected concentration and compute cellsize
    z_iface3d = d["z_iface_3d"][:,:,:,1] # assume the analysis starts from t=0 or at t_0 with given vertical grid
    z_iface3d_len = length(z_iface3d[1,1,:])
    c_force = zeros(length(d["x_center"]), length(d["y_center"]), z_iface3d_len-1)
    cellsize = similar(c_force) # the constant that translates c to NumPar
    for i in eachindex(d["x_center"])
        x = d["x_center"][i]
        for j in eachindex(d["y_center"])
            y = d["y_center"][j]
            for k=1:z_iface3d_len-1
                c_force[i,j,k] = s(x,y,(z_iface3d[i,j,k]+z_iface3d[i,j,k+1])/2,d["tstart"])# 1 unit salinity per square area per meter (computed by per unit area) is represented by 10 particles
                cellsize[i,j,k] = cosd(y) / (1/2) * (z_iface3d[i,j,k+1]-z_iface3d[i,j,k])
            end
        end
    end
    # compute expected number of particles in each cell
    NumPar_force = c_force .* cellsize  .* 10 # std (0.5*dx)*dy*1m --> 10 part
                                              # true(cos(theta)*dx)*dy*z (multiply cos and divide a halve)
    # compute true number of particles in each cell
    NumPar_real = count_particles(d["x_center"],d["y_center"],z_iface3d,p,d["isactivelist"])

    # release particles accordingly
    is_release_this_cell = (NumPar_force .- NumPar_real) .> 1e-10 
    for i in eachindex(d["x_center"])
        x = d["x_center"][i]
        for j in eachindex(d["y_center"])
            y = d["y_center"][j]
            for k=1:z_iface3d_len-1
                if is_release_this_cell[i,j,k]
                    NumPar_release = ceil(Int, NumPar_force[i,j,k] - NumPar_real[i,j,k])
                    id = d["first_inactive_index"]
                    indices = (id:id+NumPar_release-1)
                    p[1,indices] = x-d["dx"]/2 .+ d["dx"]*rand(randpool, npart_release,1) 
                    p[2,indices] = y-d["dy"]/2 .+ d["dx"]*rand(randpool, npart_release,1)
                    p[3,indices] = z_iface3d[i,j,k] .+ (z_iface3d[i,j,k+1]-z_iface3d[i,j,k])*rand(randpool, npart_release,1)
                    p[4,indices] .= 0
                    d["first_inactive_index"] += npart_release                    
                    # check if any particles are outside; if it is outside, don't activate it?
                    for l in indices
                        temp, _ = is_particle_in_domain(p[:,l])
                        if temp
                            d["isactive_particles"][l] = true
                        end
                    end

                end
            end
        end
    end    
    d["first_active_index"] = findfirst(d["isactive_particles"]) # to be commented
    return p
end
d["initialize_particles"] = initialize_particles