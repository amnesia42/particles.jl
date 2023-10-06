using DelimitedFiles, Printf, LegendrePolynomials
using Plots,Measures

z0=0.5; N = 100000000;Tend = 0.216;scheme="euler";kernel_type="Epa"
t_obs = [0.036, 0.072, 0.108, 0.144, 0.180, 0.216]
#t_obs = 3e-3:3e-3:0.216
dt_list = [3e-3,1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]
#zgrid = 0.01:0.02:0.99
zgrid = 0:0.02:1.0
Nz = length(zgrid)

TV = zeros(length(t_obs))

# plot analytical solution 
C = zeros(Nz, length(t_obs))
for j = 1:length(t_obs)
    t = t_obs[j]
    for n=0:50
        for i=1:Nz
            δC= (2*n+1)*Pl(2*zgrid[i]-1, n)*Pl(2*z0-1,n)*exp(-6*n*(n+1)*t)
            C[i,j] += δC
        end
    end
    s = sum(abs.(C[:,j] .- 1)) / Nz
    println(s)
    TV[j] = s
end
# collect the result
fname=@sprintf("./pdfoutput_kernel=%s_scheme=%s_z0=%.2f_N=%d.txt", kernel_type,scheme, z0, N)
pdfoutput = readdlm(fname,'\t',Float64,'\n')
numdt = length(dt_list)
numTobs = length(t_obs)
pdfoutput = reshape(pdfoutput, numdt, Nz, numTobs)
#pall = plot(layout=(2,ceil(Int, numdt/2)),size=(1800,900),margin=10Plots.mm)
pall = []
#pall = plot(layout=numdt,size=(1800,900),margin=10Plots.mm)
for j=1:numdt
    # compute the global error @ certain dt and t_obs
    r1 = pdfoutput[j,:,:] - C
    r2 = sum(abs, r1, dims=1)
    #r2 = r2 ./ length(r2)
    p = plot(t_obs, vec(r2), linecolor=:red, linewidth=3, label="global error", margin=10Plots.mm)
    # compute the local error at the boundary @ certain dt, t_obs
    r3 = abs.(pdfoutput[j,1,:] - C[1,:])
    plot!(t_obs, vec(r3), linecolor=:blue, linewidth=3, label="boundary local error",xlim=(0,0.3), xlabel="t", ylabel="L1error", legend=:right, dpi=300)
    # plot the percentage of contribution
    scatter!(twinx(), t_obs, vec(r3) ./ vec(r2) .* 100, linewidth=3, xlim=(0,0.3),ylim=(0,100), ylabel="percentage",label="pct.",title="dt=$(dt_list[j])",dpi=300)
    #scatter!(twinx(), t_obs, vec(r3) ./ vec(r2) .* 100, linewidth=3, label="Local error percentage", xlim=(0,0.3),ylim=(0,100), ylabel="percentage",title="dt=$(dt_list[j])", legend= :outertop, dpi=300)
    # plot the ratio
    #plot!(twinx(), t_obs, vec(r3) ./ (vec(r2) ./ Nz), markershape= :circle,linewidth=3, xlim=(0,0.3),ylim=(0,50), ylabel="ratio",label="ratio",title="dt=$(dt_list[j])",dpi=300)
    push!(pall, p)
end
pall = plot(pall..., layout=(2,ceil(Int, numdt/2)),size=(1800,900),margin=20Plots.mm,dpi=300)
figname = @sprintf("./err_boundary_kernel=%s_scheme=%s_z0=%.2f_N=%d.png",kernel_type, scheme, z0, N)
savefig(pall,figname)