# csv files with some ship tracks from AIS through the English Channel on Feb 11-15 2019 and one on
# Feb 27.
# All tracks are eastbound.

using DataFrames
using CSV
using Plots
using Dates
using Particles

# read track data
t=[]
push!(t,CSV.read("East_Bound_combined_229068000.csv",header=2))
push!(t,CSV.read("East_Bound_combined_241463000.csv",header=2))
push!(t,CSV.read("East_Bound_combined_249420000.csv",header=2))
push!(t,CSV.read("East_Bound_combined_256735000.csv",header=2))
push!(t,CSV.read("East_Bound_combined_25903000.csv",header=2))
# convert times from strings to DateTime values
for ti=1:length(t)
   times=t[ti][:,1]
   ts=[DateTime(t,"dd/mm/yyyy HH:MM:SS") for t=times]
   t[ti].Times=ts
end
#get names of the ships
shipnames=[]
for si=1:5
   push!(shipnames,t[si].Name[20])
end

# read GTSM current data
dflow_map=load_nc_info(".",r"gtsm_fine_...._map.nc");
t0=get_reftime(dflow_map)
# create interpolation functions u1(x,y,z,t) and v1(x,y,z,t) with t in seconds relative to t0
interp=load_dflow_grid(dflow_map,50,false);
u1,v1=initialize_interpolation(dflow_map,interp,t0);
# some tests
#ind=find_index(interp,-2.0,50.0)[2]
#uu=dflow_map[1].vars["ucx"];
#tu=dflow_map[1].vars["time"];
#uu[ind,1:5]
#tu[1:5]
#ff=u1(-2.0,50.0,0.0,1209600.0) #lon,lat,dummy,t with t in seconds relative to t0
# create timeseries at one point
lon=-2.0
lat=50.0
ts=(14.0*24.0*3600.0):1800.0:(18.0*24.0*3600.0) #times in seconds since t0
tt=t0+ts.*Second(1) #convert to DateTime
uu=zeros(length(tt))
vv=zeros(length(tt))
for ti=1:length(tt)
  uu[ti]=u1(lon,lat,0.0,ts[ti])
  vv[ti]=v1(lon,lat,0.0,ts[ti])
end
plot(tt,[uu,vv],label=["u east","v north"])
title!("Tidal currents from GTSM model at lon-2 lat=50")
savefig("timseries_uv_2w_50n.png")
# first interpolate to track positions and times for track 1
utrack1=zeros(length(t[1].Times))
vtrack1=zeros(length(t[1].Times))
ttrack1=zeros(length(t[1].Times))
for ti=1:length(ttrack1)
   time=(t[1].Times[ti]-t0).value/1000 #convert to seconds relative to t0
   #println("time=$(time) $(ti)")
   ttrack1[ti]=time
   utrack1[ti]=u1(t[1].Lon[ti],t[1].Lat[ti],0.0,ttrack1[ti])
   vtrack1[ti]=v1(t[1].Lon[ti],t[1].Lat[ti],0.0,ttrack1[ti])
end

#
# plot tracks on a background
#
width=1000
height=1000
Plots.default(:size,[width,height])
gebco_server=WmsServer("gebco") #gebco or emodnet-bathymetry or open-streetmap
plotbox=[-5.5,48.5,1.5,51.5] #area to plot min(Lon), min(Lat), max(Lon), max(Lat)

img=get_map(gebco_server,plotbox,width,height)

plot_image(img,plotbox)

plot!(t[1].Lon,t[1].Lat)
plot!(t[2].Lon,t[2].Lat)
plot!(t[3].Lon,t[3].Lat)
plot!(t[4].Lon,t[4].Lat)
plot!(t[5].Lon,t[5].Lat)

savefig("track_East_Bound_combined.png")

#
# plot speed as timeseries (SOG= Speed Over Ground)
#
knots2ms=0.514444 #conversion of speed as knots to meters/second
#plot(t[1].Times,t[1].SOG*knots2ms)
plot(t[1].Times,[t[1].SOG*knots2ms,utrack1,vtrack1],label=["SOG from AIS","GTSM u east","GTSM v north"])
title!("Speed and current along track1")
ylabel!("[m/s]")

savefig("timeseries_speed_1.png")
