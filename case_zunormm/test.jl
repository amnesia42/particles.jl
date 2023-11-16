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

Zarr_data = ZarrData(".","ZUNORMM_map_fullgrid_v2.zarr")
t0=get_reftime(Zarr_data)
h=initialize_interpolation(Zarr_data,"waterlevel",t0)
depth = initialize_interpolation(Zarr_data,"waterdepth",t0)
u=initialize_interpolation(Zarr_data,"x_velocity",t0)
v=initialize_interpolation(Zarr_data,"y_velocity",t0)
s=initialize_interpolation(Zarr_data,"salinity",t0)


# find the extent of the data
xy_bbox=Zarr_data.xy_grid.bbox
println("bbox=$(xy_bbox)") #bbox=[3.4773802907221345, 5.566420993890162, 51.66065856527758, 52.466444933890344]

t=Zarr_data.file.arrays["time"][:] #"seconds since 2022-04-01 00:00:00 +00:00"
println(t0) #"seconds since 2022-04-16 00:00:00 +00:00"
println( (t[2]-t[1])/3600.0 ) #2022-04-16 step=2hrs until 2022-05-01
println( (t[end]-t[2])/3600.0/24.0) #2022-04-16 duraton~=15hrs until 2022-05-01

x_left = 3.4
x_right = 5.6
y_bot = 51.6
y_top = 52.5
z_surface = 0.0
z_bottom = -40.0
t_start=0.0
t_step=2*3600.0
#t_stop=15*24*3600.0 # end of dataset at 15 days
t_stop=12*3600.0 # can be 12 hours

# This is using the interpolation functions instead of plotting the data directly, which is much faster.
# anim = @animate for t in range(t_start, stop = t_stop, step = t_step)
#     if t%1200==0 print("t=$(t) ") end
#     # grid for plotting only
#     x_points = range(x_left,stop=x_right,length=100) #only show 100 points in x
#     y_points = range(y_bot,stop=y_top,length=100)
#     z        = -0.0
#     #u_interp = [u(x,y_middle,z,t) for x in x_points, z in z_points] # x velocity
#     h_interp = [h(x,y,z,t) for x in x_points, y in y_points] # water level
#     s_interp = [s(x,y,z,t) for x in x_points, y in y_points] # salinity
#     depth_interp = [depth(x,y,z,t) for x in x_points, y in y_points] # water depth

#     l = @layout([a; b; c])
#     p1=heatmap(x_points,y_points,h_interp',xlabel="x",ylabel="y",title="waterlevel at t=$(t/3600)h",clims=(-1.5,1.5))
#     p2=heatmap(x_points,y_points,s_interp',xlabel="x",ylabel="y",title="salinity at t=$(t/3600)h",clims=(0,35))
#     plot(p1,p2,p3,layout=l,size=(800,2400))
# end
# gif(anim, "./temp/rmm_test.gif", fps = 3)

# direct plot of surface salinity at first time step - faster and nicer
x=Zarr_data.file.arrays["x_center"][:]
y=Zarr_data.file.arrays["y_center"][:]
println("length of x-coordinates=",length(x),",length of y coordinates=",length(y))
#println(Zarr_data.file.arrays["salinity"]) #salinity taken at center 
#ZArray{Int16} of size 999 x 499 x 46 x 181
#println(Zarr_data.file.arrays["z_iface_3d"])
#ZArray{Int16} of size 999 x 499 x 47 x 181

# plot waterlevel range [-1.5, 17.43]
# for tindex in range(1, stop=12, step=1)
#     waterlevel_search = 0.01 .* Array(Zarr_data.file.arrays["waterlevel"][:,:,tindex]) #999 x 499 x 181
#     waterlevel_search[waterlevel_search .>= 99.0] .= 0.0
#     @printf "minimum water level: %6.2f,\t maximum water level: %6.2f. \n" minimum(waterlevel_search) maximum(waterlevel_search)
# end

# # plot water level variation
# anim = @animate for tindex in range(1, stop=12, step=1)
#     waterlevel = 0.01 .* Array(Zarr_data.file.arrays["waterlevel"][:,:,tindex]) #999 x 499 x 181
#     waterlevel[waterlevel .>= 99.0] .= NaN
#     t = (tindex-1)*2 
#     heatmap(x,y,waterlevel',title="waterlevel at tref=$(t)h.",clims=(-1.5,3.0))
# end
# gif(anim, "./temp/waterlevel_test.gif", fps = 3)

# at the starting time
# seabed range: [-42.69, 17.43]
# sealevel range: [-0.86, 17,43] 
# waterdepth range: [0,42.82]
# waterdepth = 0.01 .* Array(Zarr_data.file.arrays["waterdepth"][:,:,1]) #999 x 499 x 181
# t = copy(waterdepth)
# t[t .>= 99.0] .= 0
# @printf "minimum waterdepth: %6.2f,\t maximum waterdepth: %6.2f. \n" minimum(t) maximum(t)
# waterlevel = 0.01 .* Array(Zarr_data.file.arrays["waterlevel"][:,:,1]) #999 x 499 x 181
# t = copy(waterlevel)
# t[t .>= 99.0] .= 0
# @printf "minimum waterlevel: %6.2f,\t maximum waterlevel: %6.2f. \n" minimum(t) maximum(t)
# r = (waterdepth - waterlevel) .* (-1)
# @printf "minimum z_bottom: %6.2f,\t maximum z_bottom: %6.2f. \n" minimum(r) maximum(r)
# waterdepth[waterdepth .>= 99.0] .= NaN
# waterlevel[waterdepth .>= 99.0] .= NaN
# r = (waterdepth - waterlevel) .* (-1)
# p = heatmap(x,y,r',title="z_{seabed}",clims=(-18,10),dpi=600)
# savefig(p, "./temp/z_seabed_distribution.png")


# # plot surface salinity
# anim = @animate for tindex in range(1, stop=12, step=1)
#     s_surface = 0.01 .* Array(Zarr_data.file.arrays["salinity"][:,:,46,tindex]) #999 x 499 x 46 x181
#     s_surface[s_surface .>= 99.0] .= NaN
#     t = (tindex-1)*2 
#     heatmap!(x,y,s_surface',title="surface salinity at tref=$(t)h.",clims=(0,40))
# end
# gif(anim, "./temp/s_surface.gif", fps = 3)

# plot bottom salinity
# anim = @animate for tindex in range(1, stop=12, step=1)
#     s_bottom = 0.01 .* Array(Zarr_data.file.arrays["salinity"][:,:,1,tindex]) #999 x 499 x 46 x181
#     s_bottom[s_bottom .>= 99.0] .= NaN
#     t = (tindex-1)*2 
#     heatmap!(x,y,s_bottom',title="bottom salinity at tref=$(t)h.",clims=(0,40))
# end
# gif(anim, "./temp/s_bottom.gif", fps = 3)
    
 
# # plot multiple z_center used in FVM computation
z_center3d = 0.01 .* Array(Zarr_data.file.arrays["z_center_3d"][:,:,:,1]) #999 x 499 x 46 x 181
z_center3d[z_center3d .>= 99.0] .= NaN

# # plot multiple z_surface used in FVM computation
z_iface3d = 0.01 .* Array(Zarr_data.file.arrays["z_iface_3d"][:,:,:,1]) #999 x 499 x 47 x 181
z_iface3d[z_iface3d .>= 99.0] .= NaN
all_surface=[]
for i=1:10
    if i==1
        p = PlotlyJS.surface(z=z_iface3d[:,:,i],x=x,y=y)
    else
        p = PlotlyJS.surface(z=z_iface3d[:,:,i],x=x,y=y,showscale=false,opacity=0.9,dpi=300)
    end
    push!(all_surface, copy(p))
end
PlotlyJS.plot(p)
savefig(p, "./temp/z_surface_distribution.png")


@printf "(x[251],y[273])=(%.2f,%.2f)" x[251] y[273]
print(z_center3d[251,273,:])
@printf "(x[500],y[273])=(%.2f,%.2f)" x[500] y[273]
print(z_center3d[500,273,:])
