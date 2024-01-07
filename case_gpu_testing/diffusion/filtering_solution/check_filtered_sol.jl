using DelimitedFiles, SpecialFunctions, LegendrePolynomials, Printf
using Plots, LaTeXStrings

kernel_type = "Epa"; scheme="euler"; z0=0.5; N=100000000

ugrid = 0:0.02:1.0
t_obs_list = [0.036, 0.072, 0.108, 0.144, 0.180, 0.216]
dt_list = [3e-3,1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]

Nz = length(ugrid)
numTobs = length(t_obs_list)
numdt = length(dt_list)

TV = zeros(length(t_obs_list))

# plot analytical solution 
C = zeros(Nz, length(t_obs_list))
for j = 1:length(t_obs_list)
    t = t_obs_list[j]
    for n=0:50
        for i=1:Nz
            δC= (2*n+1)*Pl(2*ugrid[i]-1, n)*Pl(2*z0-1,n)*exp(-6*n*(n+1)*t)
            C[i,j] += δC
        end
    end
    s = sum(abs.(C[:,j] .- 1)) / Nz
    println(s)
    TV[j] = s
end


# p_tilde of size numdt * Nz * numTobs
p_tilde = readdlm("pdfoutput_kernel=Epa_scheme=euler_z0=0.50_N=100000000.txt")
p_tilde = reshape(p_tilde, numdt, Nz, numTobs)
d= readdlm("filtered_solution_scheme=euler_z0=0.50_N=100000000_h=2.0e-05.txt", '\t')
p_filter = similar(p_tilde)
for i in 1:numdt
    p_filter[i,:,:] = reshape(d[i,:], Nz, numTobs)
end

# calculate particle-tracking solution
pall = plot(layout=(2,Int(numTobs/2)),size=(1800,900),margin=10Plots.mm)
for i=1:numTobs
    for j=1:numdt
        r =  p_tilde[j,:,i] 
        plot!(pall,r, ugrid, linewidth=3, xlim=(0,1.5),label="dt=$(dt_list[j])",dpi=300,legend=:right, legnedfontsize=16, subplot=i)
    end
    plot!(pall, subplot=i, xlabel=L"\tilde{p}(x,t)", ylabel="x",title="numerical sol at t=$(t_obs_list[i])", titlefontsize=18,xguidefontsize=18, yguidefontsize=16)
end
figname = @sprintf("./pdf_distri_kernel=%s_scheme=%s_z0=%.2f_N=%d.png",kernel_type, scheme, z0, N)
savefig(pall,figname)

# calculate filtered solution
pall = plot(layout=(2,Int(numTobs/2)),size=(2000,900),margin=10Plots.mm)
for i=1:numTobs
    for j=1:numdt
        r =  p_filter[j, :, i] 
        plot!(pall,r, ugrid, linewidth=3, marker=:circle, xlim=(0,1.5),label="dt=$(dt_list[j])",dpi=300,legend=:right, legnedfontsize=16, subplot=i)
    end
    plot!(pall, subplot=i, xlabel=L"(K\circ p)(x,t)", ylabel="x",title="filtered analytical sol at t=$(t_obs_list[i])", titlefontsize=18,xguidefontsize=18, yguidefontsize=16)
end
figname = @sprintf("./filtered_anapdf_distri_kernel=%s_scheme=%s_z0=%.2f_N=%d.png",kernel_type, scheme, z0, N)
savefig(pall,figname)

# calculate non-filtering error
pall = plot(layout=(2,Int(numTobs/2)),size=(1800,900),margin=10Plots.mm)
for i=1:numTobs
    for j=1:numdt
        r = p_filter[j,:,i] .- p_tilde[j,:,i] 
        plot!(pall,r, ugrid, linewidth=3, xlim=(-0.5,0.4),label="dt=$(dt_list[j])",dpi=300,legend=:right, legnedfontsize=16, subplot=i)
    end
    plot!(pall, subplot=i, xlabel=L"(K\circ p - E\hat{p})(x,t)", ylabel="x",title="non-filtering error at t=$(t_obs_list[i])", titlefontsize=18,xguidefontsize=18, yguidefontsize=16)
end
figname = @sprintf("./nonfiltering_err_distri_kernel=%s_scheme=%s_z0=%.2f_N=%d.png",kernel_type, scheme, z0, N)
savefig(pall,figname)

# calculate filtering error
pall = plot(layout=(2,Int(numTobs/2)),size=(1800,900),margin=10Plots.mm)
for i=1:numTobs
    for j=1:numdt
        r =  C[:,i] .- p_filter[j,:,i] 
        plot!(pall,r, ugrid, linewidth=3, marker=:circle, xlim=(-0.1,0.6),label="dt=$(dt_list[j])",dpi=300,legend=:right, legnedfontsize=16, subplot=i)
    end
    plot!(pall, subplot=i, xlabel=L"(p - K\circ p)(x,t)", ylabel="x",title="filtering error at t=$(t_obs_list[i])", titlefontsize=18,xguidefontsize=18, yguidefontsize=16)
end
figname = @sprintf("./filtering_err_distri_kernel=%s_scheme=%s_z0=%.2f_N=%d.png",kernel_type, scheme, z0, N)
savefig(pall,figname)