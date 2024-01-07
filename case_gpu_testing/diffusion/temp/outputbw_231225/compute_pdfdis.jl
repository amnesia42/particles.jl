using DelimitedFiles, Printf, LegendrePolynomials
using Plots

z0=0.5; N = 10000;Tend = 0.216;scheme="RK4";kernel_type="Epa"; h=2.0e-5
t_obs = [0.036, 0.072, 0.108, 0.144, 0.180, 0.216]
#dt_list = [3e-5]
dt_list =[3e-3,1e-3, 3e-4,1e-4,3e-5]

#dt_list = [3e-3,1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]

zgrid = 0.0:0.02:1.0
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
dirname = "case_gpu_testing/diffusion/temp/outputbw_231225"
fname=@sprintf("%s/pdfoutput_scheme=%s_z0=%.2f_N=%d_h=%.1e.txt", dirname, scheme, z0, N, h)
pdfoutput = readdlm(fname,'\t',Float64,'\n')
numdt = length(dt_list)
numTobs = length(t_obs)
pdfoutput = reshape(pdfoutput, numdt, Nz, numTobs)
pall = plot(layout=(2,Int(numTobs/2)),size=(1800,900),margin=10Plots.mm)
for i=1:numTobs
    for j=1:numdt
        plot!(pall,pdfoutput[j,:,i], linewidth=3, zgrid, label="num dt=$(dt_list[j])",dpi=300, legend=:left,legnedfontsize=16,subplot=i)
    end
    plot!(pall,C[:,i], zgrid, xlim=(0.4,1.4),label="analy",linewidth=4,legnedfontsize=14,subplot=i)
    plot!(pall, subplot=i, xlabel="L1 error", ylabel="x",title="concentration distribution at t=$(t_obs[i])", titlefontsize=16,xguidefontsize=18, yguidefontsize=18)
end
figname = @sprintf("%s/pdf_distri_scheme=%s_z0=%.2f_N=%d_h=%.1e.png",dirname, scheme, z0, N, h)
savefig(pall,figname)