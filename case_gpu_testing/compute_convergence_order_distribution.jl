using Plots, DelimitedFiles, Printf
using LegendrePolynomials

# concentration estimation parameter specification
z0=0.5; N = 100000000;Tend = 0.216;scheme="euler"
t_obs = [0.036, 0.072, 0.108, 0.144, 0.180, 0.216]
#t_obs = 3e-3:3e-3:0.216
dt_list = [3e-3,1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]
#t_obs = [0.036, 0.072, 0.108, 0.1728]

#dt_list = [3e-4,1e-4,3e-5,1e-5,3e-6]
dx = 0.02
xs = 0:dx:1.0
error_output = zeros(length(dt_list), length(t_obs))
bandwidth_output = zeros(length(dt_list), length(t_obs))
pdf_output = zeros(length(dt_list), length(xs),length(t_obs))
analpdf_times = zeros(length(xs), length(t_obs))
Nparticles = N

# compute analytical solution
for j = 1:length(t_obs)
    t = t_obs[j]
    for i=1:length(xs)
        analpdf_times[i,j] = 0.
        for n=0:50
            δC= (2*n+1)*Pl(2*xs[i]-1, n)*Pl(2*z0-1,n)*exp(-6*n*(n+1)*t)
            analpdf_times[i,j] += δC
        end
    end
end

#
# dt x t x zgrid
d = readdlm("case_gpu_testing\\diffusion\\error\\pdfoutput_kernel=Epa_scheme=euler_z0=0.50_N=100000000.txt", Float64)
d = reshape(d, length(dt_list), length(xs),length(t_obs))

res = zeros((length(dt_list), length(xs),length(t_obs)))
for i = 1:length(dt_list)
    res[i,:,:] = abs.(d[i,:,:] - analpdf_times)
end

fig = plot(layout=(2,Int(length(t_obs)/2)),size=(1900,800),margin=5Plots.mm)
Nz_spacing = 2
for j = 1:length(t_obs)
    for i in 1:Nz_spacing:26
        plot!(fig, dt_list, res[:,i,j], label="z=$(xs[i])",dpi=300, xscale=:log10, yscale=:log10, legend=:outerright, subplot=j, title="t=$(t_obs[j])s")
    end
end
savefig(fig, ".\\case_gpu_testing\\diffusion\\error\\con_order_dis\\errdis_spacing=$(Nz_spacing).png")

