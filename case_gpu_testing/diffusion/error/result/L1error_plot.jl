using DelimitedFiles, Printf
using Plots, LaTeXStrings
fname_list=["L1error_scheme=euler_z0=0.50_N=100000_Nrep=300.txt", "L1error_scheme=heun_z0=0.50_N=100000_Nrep=300.txt","L1error_scheme=m1_z0=0.50_N=100000_Nrep=305.txt","L1error_scheme=RK4_z0=0.50_N=100000_Nrep=300.txt"]
lgname_list = ["Euler", "Heun","M1","RK4"]
dt_list = [3e-3,1e-3,3e-4, 1e-4, 3e-5] #
t_obs = [0.036, 0.072, 0.108, 0.144]
numTobs = length(t_obs)
color_names = ["red", "blue", "gray", "green" ]
p2 = plot(layout=(2,Int(numTobs/2)),size=(1600,1200),margin=10Plots.mm)
for j=1:length(fname_list)
    fname = fname_list[j]
    lgname = lgname_list[j]
    temp = readdlm(fname, '\t', Float64, '\n')
    for t_index=1:numTobs
        L1errorlist = temp[t_index,:]
        plot!(p2,dt_list, L1errorlist,  xscale=:log10, yscale=:log10, minorgrid=true, label="", color=:gray, subplot=t_index)
        scatter!(dt_list, L1errorlist, markersize=6, label=lgname, legnedfontsize=15, subplot=t_index)
        plot!(p2, xlabel=L"\Delta t", ylabel=L"\frac{1}{N_z}\sum_{i=1}^{N_z}\epsilon(z_i, T_{obs})", 
        title="L1 error at T=$(t_obs[t_index])s", subplot=t_index, dpi=300,xguidefontsize=13,yguidefontsize=13)
        for i =1:length(dt_list) - 1
            tg = abs( (log10(L1errorlist[i+1]) - log10(L1errorlist[i]) ) / ( log10(dt_list[i+1]) - log10(dt_list[i]) )  )
            tg_text = @sprintf("%.3f", tg)
            if lgname=="Heun"
                #annotate!(0.5*(dt_list[i]+dt_list[i+1]), 0.3*(L1errorlist[i]+L1errorlist[i+1]), tg_text)
            else
                annotate!(0.5*(dt_list[i]+dt_list[i+1]), 0.5*(L1errorlist[i]+L1errorlist[i+1]), tg_text, subplot=t_index)
            end
        end
    end
end
savefig(p2, "./1Ddiffusion_weakerror_old.png" )