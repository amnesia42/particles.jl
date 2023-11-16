# Initialize the notebook
using Pkg
Pkg.activate("..")
#include("../src/Particles.jl")
using Particles
using Zarr
using ZipFile
using DiskArrays
using Printf
using Plots

Zarr_data = ZarrData(".","ZUNORMM_map_fullgrid_v2.zarr")
x=Zarr_data.file.arrays["x_center"][:]
y=Zarr_data.file.arrays["y_center"][:]
println("length of x-coordinates=",length(x),",length of y coordinates=",length(y))

# # plot multiple z_center used in FVM computation
z_center3d = 0.01 .* Array(Zarr_data.file.arrays["z_center_3d"][:,:,:,1]) #999 x 499 x 46 x 181
z_center3d[z_center3d .<= -99.0] .= NaN

# # plot multiple z_surface used in FVM computation
z_iface3d = 0.01 .* Array(Zarr_data.file.arrays["z_iface_3d"][:,:,:,1]) #999 x 499 x 47 x 181
z_iface3d[z_iface3d .<= -99.0] .= NaN
p = heatmap(x,y,z_iface3d[:,:,47]')
savefig(p, "./temp/ziface_surface.png")

@printf "(x[251],y[273])=(%.2f,%.2f)" x[251] y[273]
print(z_center3d[251,273,:])
@printf "(x[500],y[273])=(%.2f,%.2f)" x[500] y[273]
print(z_center3d[500,273,:])
