# Shared routines for particle modelling.
# No routines for specific data sources or particle model-equations should be included,
# just the generic stuff goes here.

using Plots
using StaticArrays
using Dates
using Printf
using LabelledArrays

debuglevel = 1

"""
   d = default_userdata()
   initialize Dict with some default configuration data
"""
function default_userdata()
    Dict(
        # general
        "dt" => 0.01,  # TODO a fixed timestep will not work in general,
        "tstart" => 0.0,
        "tend" => 1.0,
        "reftime" => DateTime(2000, 1, 1), # Jan 1st 2000,
        "coordinates" => "projected", # projected or spherical,
        "nparticles" => 10, # number of particles,
        "variables" => ("x", "y"), # recognized are x,y,lat,lon other variables are written with partial meta-data,
        "dumval" => 9999.0,
        "time_direction" => :forwards, # :forwards or :backwards,
        # plotting to screen
        "plot_maps" => false,
        "plot_maps_size" => (1200, 1000),
        "plot_maps_times" => [],
        "plot_maps_func" => plot_maps_xy,
        "plot_maps_folder" => "output",
        "plot_maps_prefix" => "map",
        # results that are kept in memmory
        "keep_particles" => false,
        "keep_particle_times" => [],
        "all_particles" => [],
        # results written to netcdf file
        "write_maps" => false,
        "write_maps_times" => [],
        "write_maps_dir" => ".",
        "write_maps_as_series" => true,
        "write_maps_filename" => "output.nc",
    )
end

# Override to quickly make
Base.getindex(nt::NamedTuple, s::AbstractString) = getindex(nt, Symbol(s))
# Base.setindex!(nt::NamedTuple, v, k::AbstractString) = setindex!(nt, v, Symbol(k))

"""
   run_simulation(d)
   Main simulation routine. The configuration is contained in Dict d,
   which may contain call-backs for plotting etc.
"""
function run_simulation(d)
    vars = d["variables"]
    npart = d["nparticles"]
    nvars = length(vars)
    p = d["particles"]
    p_all=d["all_particles"] #if requested keep particles at intermediate times
    d["particles_indices"] = []
    Plots.default(:size, d["plot_maps_size"])

    # preallocate storage for the current particle information
    p = zeros(nvars, d["nparticles_bound"])

    # show inputs
    if debuglevel > 2
        println("configuration")
        display(d)
    end

    # simulation timespan
    tstart = d["tstart"]
    tend = d["tend"]
    tref = d["reftime"]
    t = tstart
    if d["time_direction"] == :forwards
        t = tstart
    elseif d["time_direction"] == :backwards
        t = tend
    else
        throw(ArgumentError("Unsupported variable d[time_direction]"))
    end

    # initialize outputs
    target_times = Float64[]
    if d["plot_maps"]
        # init timing
        plot_maps_times = d["plot_maps_times"]
        target_times = sort(union(plot_maps_times, target_times))
        print("plotting output at t=")
        print_times(tref, plot_maps_times)
        # init plots
        Plots.default(:size, d["plot_maps_size"])
        fig1 = d["plot_maps_background"](d)
        # scatter!(fig1,p[1,:],p[2,:],legend=false)
        d["plot_maps_func"](fig1, d, p)
        # TODO possibly label=[string(t)])
        # gui(fig1)
        # check outputfolder
        plot_maps_folder = d["plot_maps_folder"]
        length(plot_maps_folder) > 0 || error("empty plot_maps_folder")
        if isdir(plot_maps_folder)
            println("Removing existing output in folder $(plot_maps_folder)")
            rm(plot_maps_folder, recursive = true)
        end
        mkdir(plot_maps_folder)
    end
    if d["write_maps"]
        write_maps_times = d["write_maps_times"]
        target_times = sort(union(write_maps_times, target_times))
        print("writing output to netcdf at t = ")
        print_times(tref, write_maps_times)
        (nc_out, ncvars) = initialize_netcdf_output(d)
    end
    if d["keep_particles"] # in memmory storage
        keep_particle_times = d["keep_particle_times"]
        target_times = sort(union(keep_particle_times, target_times))
        print("writing output to memory at t = ")
        print_times(tref, keep_particle_times)
    end

    # if the end time of the simulation is after the last output request
    # then still simulate until end times. TODO This is debatable.
    if ((length(target_times) == 0) || (target_times[end] < tend))
        push!(target_times, tend)
    end
    # remove output requests outside the simulation time-span
    if target_times[end] > tend
        temp = sort(union(target_times, tend))
        i_last = findlast(x -> x <= tend, temp)
        target_times = temp[1:i_last]
    end
    print("interrupt simulation for output at t = ")
    print_times(tref, target_times)
    println("Simulation from time $(t)s to $(tend)s since $(tref)")
    if d["time_direction"] == :forwards
        # nothing
    elseif d["time_direction"] == :backwards
        target_times = sort(target_times, rev = true)
    end


    # simulate in pieces until next output-action
    for t_stop in target_times
        p = d["particles"]
        t_abs = tref + Second(round(t))
        t_stop_abs = tref + Second(round(t_stop))
        if d["time_direction"] == :forwards
            println("t=$(t) -> $(t_stop)  : $(t_abs) -> $(t_stop_abs) : $(round(100.0 * (t_stop - tstart) / (tend - tstart), digits = 1))%")
        elseif d["time_direction"] == :backwards
            println("t=$(t) -> $(t_stop)  : $(t_abs) -> $(t_stop_abs) : $(round(100.0 * (tend - t_stop) / (tend - tstart), digits = 1))%")
        end
        t = simulate!(p, t, t_stop, d)
        println("simulate! runs for once.")
        if (d["plot_maps"]) && (t_stop in plot_maps_times)
            (debuglevel > 1) && println("plotting map output")
            Plots.default(:size, d["plot_maps_size"])
            fig1 = d["plot_maps_background"](d)
            d["plot_maps_func"](fig1, d, p)
            # sleep(1)
            title!(fig1, "time $(t_stop_abs) : t=$(t_stop)")
            # gui(fig1) #force display to screen
            prefix = d["plot_maps_prefix"]
            savefig(fig1, joinpath(d["plot_maps_folder"], @sprintf("%s_%010.2f.png", prefix, t)))
        end
        if (d["write_maps"]) && (t_stop in write_maps_times)
            write_maps_as_series = d["write_maps_as_series"]
            timei = findfirst(x -> x == t_stop, write_maps_times)
            for vari = 1:nvars
                varname = vars[vari]
                if write_maps_as_series
                    start = [timei, 1] # part x time
                    count = [1, npart]
                    NetCDF.putvar(ncvars[vari], collect(p[vari, :]'); start = start, count = count)
                else
                    start = [1, timei] # part x time
                    count = [npart, 1]
                    NetCDF.putvar(ncvars[vari], p[vari, :]; start = start, count = count)
                end
            end
        end
        if (d["keep_particles"]) && (t_stop in keep_particle_times)
            # share a reference to d["all_particles"], which is accessible by the user
            # create a collection where each element includes particle information at that time instant
            indices = findall(d["isactive_particles"])
            push!(d["particles_indices"], copy(indices))
            push!(p_all,copy(p[:,indices])) 
        end
    end

    if d["write_maps"]
        # NetCDF.close(nc_out) #close was abandoned by NetCDF
        finalize(nc_out)
    end

    # wait for user
    # if !isinteractive() #wait for user to kill final plot
    #   println("Type [enter] to finish script")
    #   readline()
    # end
end

"""
   t_next=simulate!(p,t_now,t_stop,d)

Compute multiple timesteps until t>=t_stop for all particles in matrix p.
Possibly using parameters (eg dt) from d of type userdata.
"""
function simulate!(p, t, t_stop, d)
    Δt = d["dt"]
    f! = d["f"]
    # add a flow function g that accounts for the effect of particle drift due to diffusion
    g! = d["g"]
    variables = d["variables"]
    (m, n) = size(p) # no variables x no particles
    N_activeparts = count(d["isactive_particles"])
    println("Before this updating, the number of active particle is $(N_activeparts).")
    N_releaseparts = count(p[4,:] .> 1e-10) # count the number of particles that has been released
    println("Before this updating, the number of released particles is $(N_releaseparts)")
    ∂s = @LArray zeros(length(variables)) (:x, :y, :z, :t)
    if d["time_direction"] == :forwards
        if abs(t-t_stop) < 1e-15
            # initialization
            # currently only at the boundary
            # reasonbaly the whole domain should be iterated through and assign to each location correct number of particles
            d["particles"] = d["initialize_particles"](d, p)
            d["all_particles"] = push!(d["all_particles"],copy(d["particles"]))
            if d["is_apply_cross_bc_test"] && (d["nparticles"]==1)
                d["the_virtual_trajectory"] = reshape(d["particles"], length(d["particles"]), 1)
            end
            println("Particles initialization finish!")
        end
        while (t < (t_stop - 0.25 * Δt))
            #   (debuglevel >= 2) && println("... t=$(t) < $(t_stop)")

            # add particles to or remove them from the field at appropritate times
            # adding particles is achieved by activating them and modifying its coordinates
            p = d["release_particles"](d, p, t)
            for i = d["first_active_index"]:d["first_inactive_index"]-1
            # for i = 1:n # old code iterates through all particles
                if d["isactive_particles"][i]
                    s = @view p[:, i]
                    #forward!(f!, Δt, ∂s, s, t, i, d)  # old code without drift term g
                    forward!(f!, g!, Δt, ∂s, s, t, i, d)
                end
            end
            t += Δt
        end
    elseif d["time_direction"] == :backwards
        while (t > (t_stop + 0.25 * Δt))
            #   (debuglevel >= 2) && println("... t=$(t) > $(t_stop)")
            for i = 1:n
                s = @view p[:, i]
                f!(∂s, s, t, i, d)

                # I am still not quite sure of we should use += or -=
                # I think += is the way to go, and handle the velocity difference in the f(ds,s,t,i,d) function
                s .+= ∂s .* Δt # Euler forward
                # println("   $(i) $(s)")
                #  p[:,i] = s[:]
            end
            t -= Δt
        end
    end
    return t
end

function forward!(f!, Δt, ∂s, s, t, i, d)
    f!(∂s, s, t, i, d)
    s .+= ∂s .* Δt
end

function forward!(f!, g!, Δt, ∂s, s, t, i, d)
    if !d["isactive_particles"][i]
        # check if the ith particle has been released or 
        # whether it moves out of the computation domain
        # Case 5: the particle is not yet released into the domain 
        # Case 4: the particle is already outside, no need to update it 
        # do nothing
    else
        # preallocation
        variables = d["variables"]
        δs = @LArray zeros(length(variables)) (:x, :y, :z, :t)
        ∂s_d = @LArray zeros(length(variables)) (:x, :y, :z, :t)
        ∂s_s = @LArray zeros(length(variables)) (:x, :y, :z, :t)
        s_temp = @LArray zeros(length(variables)) (:x, :y, :z, :t)

        # compute the deterministic and stochastic drift
        f!(∂s_d, s, t, i, d)
        g!(∂s_s, s, t, i, d)
        δs = ∂s_d .* Δt .+ ∂s_s .* sqrt(Δt) .* randn(d["rng"])
        s_temp = s .+ δs

        # determine if the particle is inside the boundary
        is_next_inside_domain, is_apply_reflection = d["is_particle_in_domain"](s_temp)

        if d["is_apply_cross_bc_test"] && d["isactive_particles"][i]
            # record all virtual steps in the single-particle test
            d["the_virtual_trajectory"] = hcat(d["the_virtual_trajectory"], reshape(s_temp, length(s_temp), 1))
        end
        
        if is_next_inside_domain
            # Case 1: after virtual displacement, the particle is inside the domain
            s .= s .+ δs
        elseif !d["is_cross_bc"] && is_apply_reflection
            # Case 2: after virtual displacement, the particle is outside the domain, so halve dt is needed
            # prevent the recursion algorithm with too many depths
            if Δt < d["dt"]/2^d["recursion_depth"]
            # Case 3: the dt is too small to be halved, so let the particle cross the boundary
                d["isactive_particles"][i] = false
                return nothing
            end             
            forward!(f!, g!, Δt/2, ∂s, s, t, i, d)
            forward!(f!, g!, Δt/2, ∂s, s, t, i, d)            

        else
            # case 6: the particles leaves at the Neumann boundary
            # the particle will stay outside the domain
            s .= s .+ δs
            d["isactive_particles"][i] = false
        end
    end
end


"""
print_times(reftime,times)

Print an array of relative times in compact readable format
"""
function print_times(reftime, times)
    println(IOContext(stdout, :limit => true), times)
end

"""
   i1 = index(2,[1,2,3,4])
   i2 = index("bob",["alex","bob","charlie"])
   Find first occurrence of a variable in an array.
   Returns nothing if the value is not found.
"""
function index(var, vars)
    return indexin([var], vars)[1]
end

"""
   f=plot_maps_xy(fig,d,p)
   Plot particles as dots in xy-plane.
"""
function plot_maps_xy(fig, d, p)
    if "x" in d["variables"]
        x_index = index("x", d["variables"])
        y_index = index("y", d["variables"])
    elseif "lon" in d["variables"]
        x_index = index("lon", d["variables"])
        y_index = index("lat", d["variables"])
    else
        error("plot_maps_xy: no spatial variables x,y or lat,lon found")
    end
    scatter!(fig, p[x_index, :], p[y_index, :], markercolor = :black, legend = false)
end

"""
   f=plot_maps_xz(fig,d,p)
   Plot particles as dots in xy-plane.
"""
function plot_maps_xz(fig, d, p)
    x_index = index("x", d["variables"])
    z_index = index("z", d["variables"])
    scatter!(fig, p[x_index, :], p[z_index, :], markercolor = :black, legend = false)
end

function initialize_netcdf_output(d)
    # write as series (loc,time) or maps (times,locs)
    write_maps_as_series = true # default
    if haskey(d, "write_maps_as_series")
        write_maps_as_series = d["write_maps_as_series"]
    end
    # file
    filedir = d["write_maps_dir"]
    filename = d["write_maps_filename"]
    fullfile = joinpath(filedir, filename)
    if !isdir(filedir)
        println("Directory for output does not exist $(filedir). Creating it.")
        mkpath(filedir)
    end
    if isfile(fullfile)
        println("Output file exists. Removing file $(fullfile)")
        rm(fullfile)
    end
    # time
    t0 = d["reftime"]
    # time_atts = Dict("units" => "seconds since $(Dates.format(t0,"yyyy-mm-dd HH:MM:SS"))",
    #                 "standard_name" => "time", "long_name" => "time",
    #                 "comment" => "unspecified time zone", "calendar" => "gregorian" )
    time_atts = Dict("units" => "seconds since $(Dates.format(t0, "yyyy-mm-dd HH:MM:SS"))",
        "standard_name" => "time", "long_name" => "Time",
        "comment" => "unspecified time zone", "calendar" => "gregorian")
    map_times = collect(d["write_maps_times"])
    time_dim = NcDim("time", map_times, time_atts)
    # particle
    npart = d["nparticles"]
    vars = d["variables"]
    nvars = length(vars)
    part_atts = Dict("long_name" => "particle id", "units" => "1", "cf_role" => "trajectory_id", "missing_value" => 9999)
    part_dim = NcDim("particles", collect(1:1:npart), part_atts)

    # global attributes
    gatts = Dict("title" => "Output of particle model", "Conventions" => "CF-1.6", "featureType" => "trajectory")

    myvars = []
    dumval = d["dumval"]
    for vari = 1:nvars
        varname = vars[vari]
        varatts = Dict("long_name" => varname, "missing_value" => Float64(dumval))
        if varname == "x"
            varatts["long_name"] = "x-coordinate"
            varatts["standard_name"] = "projection_x_coordinate"
            varatts["units"] = "m"
        elseif varname == "y"
            varatts["long_name"] = "y-coordinate"
            varatts["standard_name"] = "projection_y_coordinate"
            varatts["units"] = "m"
        elseif varname == "z"
            varatts["long_name"] = "z-coordinate"
            varatts["units"] = "m"
        elseif varname == "lon"
            varatts["long_name"] = "Longitude"
            varatts["standard_name"] = "longitude"
            varatts["units"] = "degrees_east"
        elseif varname == "lat"
            varatts["long_name"] = "Latitude"
            varatts["standard_name"] = "latitude"
            varatts["units"] = "degrees_north"
        elseif varname == "age"
            varatts["long_name"] = "Age of particles"
            varatts["units"] = "s"
            varatts["coordinates"] = "time lat lon"
        else
            varatts["coordinates"] = "time lat lon"
        end
        if write_maps_as_series
            myvar = NcVar(varname, [time_dim, part_dim], atts = varatts, t = Float64)
        else
            myvar = NcVar(varname, [part_dim, time_dim], atts = varatts, t = Float64)
        end

        push!(myvars, myvar)
    end

    nc = NetCDF.create(fullfile, NcVar[myvars...], gatts = gatts, mode = NC_NETCDF4)

    p = d["particles"] # var x part
    for vari = 1:nvars
        start = [1, 1] # part x time
        if write_maps_as_series
            count = [1, npart]
            NetCDF.putvar(myvars[vari], collect(p[vari, :]'); start = start, count = count)
        else
            count = [npart, 1]
            NetCDF.putvar(myvars[vari], p[vari, :]; start = start, count = count)
        end
    end

    return (nc, myvars)
    # NetCDF.close(nc)
end
