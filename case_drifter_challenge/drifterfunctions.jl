using DelimitedFiles
using Dates

############## Put info drifters in matrix: time, longitude, latitude ##############
function drifterdata(number)
   # joinpath should be os independent, no special iswindows or \\ needed
   if Sys.iswindows()
      rawdata = readdlm("data\\drifters\\North_Sea_Drifter$(number)_GPS.tab")
   else
      rawdata = readdlm(joinpath(@__DIR__, "data/drifters/North_Sea_Drifter$(number)_GPS.tab"))                       # Get data drifter files
   end
   indexdata = findfirst(x -> x == "Longitude", rawdata)                           # Find "Longitude"
   strippeddrifter = rawdata[(indexdata[1] + 1):end,1:indexdata[2]]                # All rows after "Longitude" and all columns from 1 to the column of "Longitude" should be stored in new matrix: relevant data of the drifter

   if number <= 4 || number == 7                                                 # Because of gap in SWAN DCSM files on 10 October, drifters 5 and 6 should be started after that gap
      drifternumber = strippeddrifter
   elseif number == 5
      startdrifter = DateTime(2017, 10, 10, 19)
      startindex = findfirst(x -> x == startdrifter, floor.(DateTime.(strippeddrifter[:,1]), Dates.Hour))  # Find first time after gap
      drifternumber = strippeddrifter[startindex[1]:end,1:3]                     # Take only data after gap
   else
      startdrifter = DateTime(2017, 10, 10, 18)
      startindex = findfirst(x -> x == startdrifter, floor.(DateTime.(strippeddrifter[:,1]), Dates.Hour))  # Find first time after gap
      drifternumber = strippeddrifter[startindex[1]:end,1:3]                     # Take only data after gap
   end

   return drifternumber
end

############## Retrieve reftime, starttime and endtime of the drifter ##############
function driftertimes(drifter)
    startdrifter = DateTime(drifter[1,1])                                        # First time drifter is in water
    reftime = DateTime(Date(startdrifter))                                       # Day of which drifter is first time in water, reference time
    Hr = Dates.value(Hour(startdrifter))                                         # Take hours
    Min = Dates.value(Minute(startdrifter))                                      # Take Minutes
    Sec = Dates.value(Second(startdrifter))                                      # Take Seconds

    starttime = Hr * 3600 + Min * 60 + Sec                                           # Convert everything to hours

    enddrifter = DateTime(drifter[end,1])                                        # Time at which drifter washed ashore
    endtime = ((Dates.value(enddrifter - reftime)) / 1000)                    # Convert endtime to seconds

    return reftime, starttime, endtime
end

############## Convert the magnitude to x and y component ##############
function decompose(dir, magnitude)
   deg2rad = pi / 180.0                                                            # Converts degrees to radians

   if dir >= 270.0 && dir < 360.0
      dir = 360.0 - dir
      y_comp = -1 * cos(dir * deg2rad) * magnitude
      x_comp = 1 * cos((90.0 - dir) * deg2rad) * magnitude
   elseif dir >= 180.0 && dir < 270.0
      dir = dir - 180.0
      y_comp = 1 * cos(dir * deg2rad) * magnitude
      x_comp = 1 * cos((90.0 - dir) * deg2rad) * magnitude
   elseif dir >= 90.0 && dir < 180.0
      dir = 180.0 - dir
      y_comp = 1 * cos(dir * deg2rad) * magnitude
      x_comp = -1 * cos((90.0 - dir) * deg2rad) * magnitude
   else
      dir = dir
      y_comp = -1 * cos(dir * deg2rad) * magnitude
      x_comp = -1 * cos((90.0 - dir) * deg2rad) * magnitude
   end
   return x_comp, y_comp
end


############## Calculates JONSWAP Stokes drift ##############
function uv_sJ(Hs, Tm, direction)
   if Hs == Tm == direction == 0.0
      return 0.0, 0.0
   end
   g = 9.81                                                                      # Gravitational acceleration
   Tp = Tm * (1 / 0.95)                                                              # Conversion from mean wave period to peak period
   U_sJON = 3.18 * pi^3 / g * Hs^2 / Tm^3 * 0.4                                        # Calculating the Stokes velocity using JONSWAP.

   u_sJON, v_sJON = decompose(direction, U_sJON)                                  # Decompose Stokes velocity in x and y component

   if isnan(u_sJON) || isnan(v_sJON) || isinf(u_sJON) || isinf(v_sJON)
      u_sJON = 0.0
      v_sJON = 0.0
   end

   return u_sJON, v_sJON
end

############## Inclusion of wind drag and calculating final particle velocity ##############
function water_stokes_wind(u_a, v_a, u_w, v_w, u_s, v_s)
   	# Typical values for air and water in the North Sea at 15 degrees celcius and 35g/kg salinity
   	ρ_w = 1026                                                                 # Water density
   	ρ_a = 1.217                                                                # Air density

   	A_w_x = 0.11946                                                            # Cross-sectional area of drifter below water surface
   	A_w_y = 0.11946
   	A_a_x = 0.000623                                                           # Cross-sectional area of drifter above water surface
   	A_a_y = 0.000623

   	Cd_w = 0.63                                                                # Drag-coefficient water
   	Cd_a = Cd_w / 1.3350745019253307                                             # Drag-coefficient air

   	# Constants k of water and air
   	k_w_x = sqrt(Cd_w * ρ_w * A_w_x)
   	k_w_y = sqrt(Cd_w * ρ_w * A_w_y)

   	k_a_x = sqrt(Cd_a * ρ_a * A_a_x)
   	k_a_y = sqrt(Cd_a * ρ_a * A_a_y)

   	# Wind velocity at protruding height near sea surface
   	z_0 = 0.0002                                                               # Roughness length

   	u_a_s = u_a * (log(0.0453067 / z_0) / log(10 / z_0))				                     # Wind at protruding height near sea surface
   	v_a_s = v_a * (log(0.0453067 / z_0) / log(10 / z_0))

      uws = u_w + u_s
      vws = v_w + v_s

		u_p = (k_w_x * uws + k_a_x * u_a_s) / (k_w_x + k_a_x)
		v_p = (k_w_y * vws + k_a_y * v_a_s) / (k_w_y + k_a_y)


   	return u_p, v_p
end

############## For plotting the drifter marker ##############
function track_of_drifter!(ds, s, t, t0, t_step, drifter)
   t_sec = Second(round(t + t_step))                                               # Telling julia that the number is in seconds
   t_actual = t0 + t_sec                                                         # Convert to DateTime

   ms = t_actual .- DateTime.(drifter[:,1])
   if Dates.value(ms[1]) >= 0
     indextime = findfirst(x -> x == minimum(abs.(ms)), abs.(ms))                   # Finding time in the dataset that is closest to the present simulation time
     if indextime != nothing                                                     # If simulation time is present in the drifter data file. If the time is not present, nothing changes. So the values of the previous time step will be used.
        s[1] = drifter[indextime,3]                                            # Particle gets same coordinates as the drifter at that time
        s[2] = drifter[indextime,2]
     end
   else
     println("Take starttime of run after starttime of the drifter")
   end
   ds[1] = 0;
   ds[2] = 0;

   ds[3] = dage = 1.0

   return ds, s
end

"""
   (ν, dν/dx, dν/dy) = estimate_viscosity_smag(interp,x,y,t,u,v,default_value=1.0; vb=0.1,Sc=0.2,calc_derivatives=true)

   Estimates the Eddy viscosity and its spatial derivates at a given point (x,y,t)
   By default, it uses a background viscosity of 0.1 and a Smagorinsky coefficient of 0.1

   This function is very inefficient as it has to interpolate the water velocities at different positions, to estimate the derivatives
   If used more often, this should be sped up, e.g. use neighbouring cells. Now, using this function for each particle i, the runtime
   is slowed down by approximately for times, as there are 4 times as much interpolation request f(x,y,z,t)
"""
function estimate_viscosity_smag(interp, x, y, t, u, v, default_value=1.0; vb=0.1,Sc=0.1,calc_derivatives=true)
   ind = find_index(interp, x, y)
   if ind[1] == -1 || ind[2] == -1
      return (default_value, 0.0, 0.0)
   end
   corners = interp.grids[ind[1]].edges[ind[2],1:interp.grids[ind[1]].nnodes[ind[2]]]
   xnode = interp.grids[ind[1]].xnodes[corners]
   ynode = interp.grids[ind[1]].ynodes[corners]
   dlon = maximum(xnode) - minimum(xnode)
   dlat = maximum(ynode) - minimum(ynode)
   z = 0.0
   u_cell = u(x, y, z, t)
   v_cell = v(x, y, z, t)
   if u_cell > 0
      u_dx = u(x + dlon, y, z, t) - u_cell
      v_dx = v(x + dlon, y, z, t) - v_cell
   else
      u_dx = u_cell - u(x - dlon, y, z, t)
      v_dx = v_cell - v(x - dlon, y, z, t)
   end
   if v_cell > 0
      u_dy = u(x, y + dlat, z, t) - u_cell
      v_dy = v(x, y + dlat, z, t) - v_cell
   else
      u_dy = u_cell - u(x, y - dlat, z, t)
      v_dy = v_cell - v(x, y - dlat, z, t)
   end
   R = 6371.0e3                                                                  # Mean radius of earth from wikipedia
   deg2rad = pi / 180.0                                                            # Converts degrees to radians
   rad2deg = 180.0 / pi
   dx = deg2rad * R * cos(ynode[1] * deg2rad) * dlon # in meters
   dy = deg2rad * R * dlat                    # in meters
   viscosity = vb + (Sc * sqrt(dx * dy))^2 * sqrt(2 * (u_dx / dx)^2 + (u_dy / dy + v_dx / dx)^2 + 2 * (v_dy / dy)^2)
   if calc_derivatives
      if u_cell > 0
         viscosity_dx = estimate_viscosity_smag(interp, x + dlon, y, t, u, v, default_value; vb=vb,Sc=Sc,calc_derivatives=false)[1] - viscosity
      else
         viscosity_dx = viscosity - estimate_viscosity_smag(interp, x - dlon, y, t, u, v, default_value; vb=vb,Sc=Sc,calc_derivatives=false)[1]
      end
      if v_cell > 0
         viscosity_dy = estimate_viscosity_smag(interp, x, y + dlat, t, u, v, default_value; vb=vb,Sc=Sc,calc_derivatives=false)[1] - viscosity
      else
         viscosity_dy = viscosity - estimate_viscosity_smag(interp, x, y - dlat, t, u, v, default_value; vb=vb,Sc=Sc,calc_derivatives=false)[1]
      end
      return (viscosity, viscosity_dx / dx, viscosity_dy / dy)
   else
      return (viscosity, 0.0, 0.0)
   end
end
nothing
