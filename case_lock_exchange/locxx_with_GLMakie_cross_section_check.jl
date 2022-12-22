# Load required libraries
using Pkg
Pkg.instantiate()
Pkg.activate("..")
using Particles
using GLMakie
GLMakie.activate!()
using Random
using NetCDF
using Printf


randpool = MersenneTwister(0) #Generate same random numbers every time

# load useful default settings 
d=default_userdata()
n=10 #number of particles

# general setting
d["coordinates"] = "projected"
d["nparticles"]= 10 #number of particles
d["time_direction"] = :forwards # :forwards or :backwards
d["variables"] = ["x","y","z","age"] # if not stated otherwise explicitly

# load data from nc file 
dflow_map = load_nc_info(@__DIR__, r"locxx_map.nc")
# construct a grid for interpolation
# 2nd argument is number of grid points per dimension
# 3rd argument to supress spherical coordinates
interp = load_dflow_grid(dflow_map, 50, false); 

# set simulation time
# read the reftime from the file
d["reftime"] = get_reftime(dflow_map)
t_ref = get_reftime(dflow_map)
rel_times = get_times(dflow_map, t_ref)
d["dt"] = rel_times[2] - rel_times[1]
d["tstart"] = rel_times[begin]
d["tend"] = rel_times[end]

# my way to load the data
function load_nc_dims(ncFile_name)
    nc = NetCDF.open(ncFile_name, readdimvar=true)
    dims = Dict()
    for k in keys(nc.dim)
        display(k)
        dims[k] = Int(nc.dim[k].dimlen)
    end
    return dims
end
ncFile_name = "locxx_map.nc"
nc = NetCDF.open(ncFile_name, readdimvar=true)
# return a dict that stores the name and length of different dimensions using user-defined functions
nc_dims = load_nc_dims("locxx_map.nc") 
# particle settings 
xnodes = nc.vars["mesh2d_node_x"]
ynodes = nc.vars["mesh2d_node_y"]
znodes = nc.vars["mesh2d_interface_z"]
xmin = minimum(xnodes); xmax = maximum(xnodes);
ymin = minimum(ynodes); ymax = maximum(ynodes);
zmin = minimum(znodes); zmax = maximum(znodes);
p0 = zeros(length(d["variables"]), d["nparticles"]) # the location of particle at this instant
p0[1, :] = 0.5*(xmin+xmax)*(1 .+ rand(randpool, d["nparticles"])) 
p0[2, :] .= 0.5*(ymin+ymax)
p0[3, :] = 0.5*(zmin+zmax)*(0 .+ rand(randpool, d["nparticles"]))
d["particles"] = copy(p0)

# record grid range in each dimension
d["bbox"] = (xmin, ymin, zmin, xmax, ymax, zmax)

# plot map times 
d["plot_maps_times"] = rel_times                        # Time at which plot should be made
d["plot_maps"] = false
d["plot_maps_folder"] = "maps"
d["plot_maps_size"]=(2000,600)

# plot particle times
d["keep_particle_times"] = rel_times                        # Time at which plot should be made
d["keep_particles"] = true

# prepare the interpolator for velocity field
# no need for the 3rd dimension 
u = initialize_interpolation(dflow_map, interp, "mesh2d_ucx", d["reftime"], 0.0, d["time_direction"]);
v = initialize_interpolation(dflow_map, interp, "mesh2d_ucy", d["reftime"], 0.0, d["time_direction"]);
w = initialize_interpolation(dflow_map, interp, "mesh2d_ucz", d["reftime"], 0.0, d["time_direction"]);

# explicitly prepare the grid for ploting
scale = 10
x_resolution=31 #[0,300]
y_resolution=31   #[0,1]
x = collect(range(xmin, xmax, length=x_resolution))
y = collect(range(ymin, ymax, length=y_resolution))
xs=kron(x,ones(length(y)))
ys=kron(ones(length(x)),y)
u_plot = u.(xs, ys, 0, 100)
v_plot = v.(xs, ys, 0, 100)
strength = vec(sqrt.(u_plot.^2 .+ v_plot.^2))
f = Figure(resolution = d["plot_maps_size"])
ax = Axis(f[1, 1])
arrows!(ax, xs, ys, u_plot, v_plot, arrowsize=5, lengthscale=10,
        arrowcolor = strength, linecolor=strength)
save("bkimg_test.png", f) 

###### Velocity function for the particles ######
function f1!(ds, s, t, i, d)
    x, y, z, age = s
    z = 0.0
    up = 0
    vp = 0
    dt = d["dt"]
    uw = u(x, y, z, t) # the second w is a subscript
    vw = v(x, y, z, t)
    ww = w(x, y, z, t)
 
    # The model with :
    # 1: Only flow velocities
    up = uw
    vp = vw
    wp = ww
    # original computation for computing the particle velocity
    #(K, Kdx, Kdy) = estimate_viscosity_smag(interp, x, y, t, u, v)
    #if !(uw == vw == ua == va == 0.0)
    #   # https://doi.org/10.1016/j.ocemod.2017.11.008 eq. 27
    #   up += Kdy + randn() * sqrt(2 * K * dt) / dt
    #   vp += Kdx + randn() * sqrt(2 * K * dt) / dt
    #end
    
    epsx = 0#1
    epsz = 0#1e-5
    up += randn()*sqrt(2*epsx*dt)/dt
    vp += randn()*sqrt(2*epsz*dt)/dt
    wp += 0
 
    ds.x = up
    ds.y = vp
    ds.z = wp
    ds.t = 1.0
 
 
    if d["time_direction"] == :backwards
       up *= -1
       vp *= -1
       wp *= -1
    end
 end
 d["f"]=f1!

 # use streamfunction as background for plotting
function plot_vfield(d, t, plane; x_resolution=41, y_resolution=21, z_resolution=21, y_typical = 0.5)
    # x_resolution=41 #[0,300]
    # y_resolution=5   #[0,1]
    # z_resolution=21   #[-10,0]
    x_min, y_min, z_min, x_max, y_max, z_max = d["bbox"];
    if plane == "xy"
        x = collect(range(x_min, x_max, length=x_resolution))
        y = collect(range(y_min, y_max, length=y_resolution))
        xs=kron(x,ones(length(y)))
        ys=kron(ones(length(x)),y)
        z_typical = 0.5 * (z_min + z_max)
        u_plot = u.(xs, ys, z_typical, t)
        v_plot = v.(xs, ys, z_typical, t)
        strength = vec(sqrt.(u_plot.^2 .+ v_plot.^2))
        f = Figure(resolution=d["plot_maps_size"])
        ax = Axis(f[1,1])
        arrows!(xs, ys, u_plot, v_plot, arrowsize=5, lengthscale=10,
                arrowcolor = strength, linecolor=strength)
    elseif plane == "xz"
        x = collect(range(x_min, x_max, length=x_resolution))
        z = collect(range(z_min, z_max, length=z_resolution))
        xs=kron(x,ones(length(z)))
        zs=kron(ones(length(x)),z)
        #y_typical = 0.25 * (y_min + y_max)
        u_plot = u.(xs, y_typical, zs, t)
        w_plot = w.(xs, y_typical, zs, t)
        strength = vec(sqrt.(u_plot.^2 .+ w_plot.^2))
        f = Figure(resolution=d["plot_maps_size"])
        ax = Axis(f[1,1])
        arrows!(xs, zs, u_plot, w_plot, arrowsize=5, lengthscale=10,
                arrowcolor = strength, linecolor=strength)
    elseif plane == "yz"
        y = collect(range(y_min, y_max, length=y_resolution))
        z = collect(range(z_min, z_max, length=z_resolution))
        ys=kron(x,ones(length(y)))
        zs=kron(ones(length(x)),y)
        x_typical = 0.5 * (x_min + x_max)
        v_plot = v.(x_typical, ys, zs, t)
        w_plot = w.(x_typical, ys, zs, t)
        strength = vec(sqrt.(v_plot.^2 .+ w_plot.^2))
        f = Figure(resolution=d["plot_maps_size"])
        ax = Axis(f[1,1])
        arrows!(ys, zs, v_plot, w_plot, arrowsize=5, lengthscale=10,
                arrowcolor = strength, linecolor=strength)
    end
    return f,ax
end
d["plot_vfield"] = plot_vfield 

function plot_particles(ax, d, p, plane)
    if plane == "xy"
        index_dir1 = index("x", d["variables"])
        index_dir2 = index("y", d["variables"])
    elseif plane == "xz"
        index_dir1 = index("x", d["variables"])
        index_dir2 = index("z", d["variables"])    
    elseif plane == "yz"
        index_dir1 = index("y", d["variables"])
        index_dir2 = index("z", d["variables"])    
    end
    scatter!(ax, p[index_dir1, :], p[index_dir2, :], color = :black)
end

# screenshot of particles
function plot_screenshot(d; dirname="", fnametag="", y_typical=0.5, plane="xz")
    t_plot = d["keep_particle_times"];
    p = d["all_particles"]; # because of "the use of keep_particles"
    for i in 1:length(t_plot)
        fname = dirname*"\\"*plane*"_screenshot_time"*"$(@sprintf("%.2f", t_plot[i]))"*".png"
        fig, ax = d["plot_vfield"](d, t_plot[i], plane; y_typical=y_typical)
        #ax = plot_particles(ax, d, p[i], plane)
        save(fname, fig)  
    end
end

open("stdout.txt","w") do io
    redirect_stdout(io) do
        for y_typical in [0.1, 0.2, 0.3, 0.4, 0.6, 0.8]
            println("y_typical=$y_typical")
            dirname = "vortex_search\\$y_typical y"
            if ~isdir(dirname)
                mkdir(dirname)
            end
            dirname = joinpath(@__DIR__, dirname)
            @time run_simulation(d)
            plot_screenshot(d; dirname=dirname, y_typical=y_typical)
        end
    end
end

