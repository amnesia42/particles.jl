using DelimitedFiles, Printf, LegendrePolynomials
using Plots

z0=0.5; N = 10000000;Tend = 0.216;scheme="euler";kernel_type="Epa"
t_obs = [0.036, 0.072, 0.108, 0.144, 0.180, 0.216]
#t_obs = 3e-3:3e-3:0.216
dt_list = [3e-3,1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]
dz = 0.02
zgrid = 0:dz:1.0
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
pall = plot(layout=(2,Int(numTobs/2)),size=(1800,900),margin=10Plots.mm)
for i=1:numTobs
    #plot!(pall,title="L1 absolute error distribution at t=$(t_obs[i])",subplot=i)
    for j=1:numdt
        plot!(pall,pdfoutput[j,:,i], linewidth=3, zgrid, xlabel="concentration", ylabel="x",label="num dt=$(dt_list[j])",dpi=300, legend=:left,subplot=i)
    end
    plot!(pall,C[:,i], zgrid, xlim=(0.4,1.4),label="analy",linewidth=4,subplot=i)
end
figname = @sprintf("./pdf_distri_kernel=%s_scheme=%s_z0=%.2f_N=%d.png",kernel_type, scheme, z0, N)
savefig(pall,figname)