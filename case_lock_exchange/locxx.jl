# Load required libraries
using Pkg
Pkg.instantiate()
Pkg.activate("..")
using Plots
using Particles
using GLMakie
using Random

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

# particle settings
xnodes = dflow_map[1].vars["mesh2d_node_x"]
ynodes = dflow_map[1].vars["mesh2d_node_y"]
xmin = minimum(xnodes); xmax = maximum(xnodes);
ymin = minimum(ynodes); ymax = maximum(ynodes);
p0 = zeros(length(d["variables"]), d["nparticles"])
p0[1, :] = 0.5*(xmin+xmax)*(1 .+ rand(randpool, d["nparticles"])) 
p0[2, :] = 0.5*(ymin+ymax)*(1 .+ rand(randpool, d["nparticles"]))
d["particles"] = copy(p0)

# record grid range in each dimension
d["bbox"] = (xmin, ymin, xmax, ymax)

# plot map times 
d["plot_maps_times"] = rel_times                        # Time at which plot should be made
d["plot_maps"] = true
d["plot_maps_folder"] = "maps"
d["plot_maps_size"]=(2000,600)

# plot particle times
d["keep_particle_times"] = rel_times                        # Time at which plot should be made
d["keep_particles"] = true

# prepare the interpolator for velocity field
# no need for the 3rd dimension 
u = initialize_interpolation(dflow_map, interp, "mesh2d_ucx", d["reftime"], 0.0, d["time_direction"]);
v = initialize_interpolation(dflow_map, interp, "mesh2d_ucy", d["reftime"], 0.0, d["time_direction"]);
display(typeof(u))

GLMakie.activate!()
# explicitly prepare the grid for ploting
scale = 10
x_resolution=31 #[0,300]
y_resolution=31   #[0,1]
x = collect(range(minimum(dflow_map[1]["mesh2d_node_x"][:]), maximum(dflow_map[1]["mesh2d_node_x"][:]), length=x_resolution))
y = collect(range(minimum(dflow_map[1]["mesh2d_node_y"][:]), maximum(dflow_map[1]["mesh2d_node_y"][:]), length=y_resolution))
xs=kron(x,ones(length(y)))
ys=kron(ones(length(x)),y)

Plots.default(:size, d["plot_maps_size"])
u_plot = zeros(length(x)*length(y))
v_plot = zeros(length(x)*length(y))
u_plot = u.(xs, ys, 0, 100)
v_plot = v.(xs, ys, 0, 100)
strength = vec(sqrt.(u_plot.^2 .+ v_plot.^2))
f = Figure(resolution = d["plot_maps_size"])
ax = Axis(f[1, 1])
arrows!(xs, ys, u_plot, v_plot, arrowsize=5, lengthscale=10,
        arrowcolor = strength, linecolor=strength)
save("bkimg.png", f) 

###### Velocity function for the particles ######
function f1!(ds, s, t, i, d)
    x, y, z, age = s
    z = 0.0
    up = 0
    vp = 0
    dt = d["dt"]
    uw = u(x, y, z, t)
    vw = v(x, y, z, t)
 
    # The model with :
    # 1: Only flow velocities
    up = uw
    vp = vw
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
 
    ds.x = up
    ds.y = vp
    ds.z = 0.0
    ds.t = 1.0
 
 
    if d["time_direction"] == :backwards
       up *= -1
       vp *= -1
    end
 end
 d["f"]=f1!

 # use streamfunction as background for plotting
function plot_vfield(d,t)
    GLMakie.activate!()
    x_resolution=41 #[0,300]
    y_resolution=21   #[0,1]
    x_min, y_min, x_max, y_max = d["bbox"]
    x = collect(range(x_min, x_max, length=x_resolution))
    y = collect(range(y_min, y_max, length=y_resolution))
    xs=kron(x,ones(length(y)))
    ys=kron(ones(length(x)),y)
    Plots.default(:size, d["plot_maps_size"])
    u_plot = u.(xs, ys, 0, t)
    v_plot = v.(xs, ys, 0, t)
    strength = vec(sqrt.(u_plot.^2 .+ v_plot.^2))
    f = Figure()
    ax = Axis(f[1,1])
    arrows!(xs, ys, us, vs, arrowsize=5, lengthscale=10,
            arrowcolor = strength, linecolor=strength)
    return(f)
 end
 d["plot_vfield"] = plot_vfield 

 # screenshot of particles
function plot_screenshot(d; dirname="", fnametag="")
    GLMakie.activate!()
    t_plot=d["keep_particle_times"];
    p=d["all_particles"]
    for i in 1:length(t_plot)
        fig = d["plot_vfield"](d, t_plot[i])
        fname = dirname*"\\screenshot_time"*"$(t_plot[i])"*".png"
        d["plot_maps_func"](fig, d, p[i])
        save(fname, fig)
    end
end

dirname = "test0"
if ~isdir(dirname)
    mkdir(dirname)
end
dirname = joinpath(@__DIR__, dirname)

open("stdout.txt","w") do io
    redirect_stdout(io) do
        @time run_simulation(d)
        plot_screenshot(d; dirname=dirname)
    end
end