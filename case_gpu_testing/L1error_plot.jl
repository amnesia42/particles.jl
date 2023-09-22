using DelimitedFiles
using Plots, LaTeXStrings
d = readdlm("./diffusion/error/result/L1error_scheme=euler_z0=0.50_N=100000_Nrep=300.txt", '\t', Float64, '\n')
L1errorlist = d[3,:] # take the results at T=0.108
dt_list = [3e-3,1e-3,3e-4, 1e-4, 3e-5] #
t_obs = [0.036, 0.072, 0.108, 0.144]

p2 = plot(dt_list, L1errorlist,  xscale=:log10, yscale=:log10, minorgrid=true, label="")
scatter!(dt_list, L1errorlist, label="L1 error")
for i =1:length(dt_list) - 1
    tg = abs( (log10(L1errorlist[i+1]) - log10(L1errorlist[i]) ) / ( log10(dt_list[i+1]) - log10(dt_list[i]) )  )
    tg_text = @sprintf("%.3f", tg)
    annotate!(0.5*(dt_list[i]+dt_list[i+1]), 0.5*(L1errorlist[i]+L1errorlist[i+1]), tg_text)
end
plot!(p2, xlabel=L"\Delta t", ylabel=L"E[TV[c_{num}(x,T)-c_{anal}(x,T)]]", title="L1 error at T=$(0.108)", dpi=300)
savefig(p2, "./Euler_1Ddiffusion_weakerror.png" )