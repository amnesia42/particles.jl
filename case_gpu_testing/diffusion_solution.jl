include("ConcentrationCalculationLibrary.jl")
using LegendrePolynomials
using Plots
z0=0.5
Nz = 101
t_list = [0.036 0.072 0.108 0.144 0.180 0.216]
#t_list = [0.003]
zgrid = range(0,1,Nz)
TV = zeros(length(t_list))

# plot analytical solution 
C = ones(Nz, length(t_list))
for j = 1:length(t_list)
    t = t_list[j]
    for n=1:50
        for i=1:Nz
            δC= (2*n+1)*Pl(2*zgrid[i]-1, n)*Pl(2*z0-1,n)*exp(-6*n*(n+1)*t)
            C[i,j] += δC
        #     if i==51
        #         println("n={$n}, C(0.5)=$(δC)")
        #      end
        #     if i==101
        #        println("n={$n}, C(1.0)=$(δC)")
        #     end
        end
    end
    s = sum(abs.(C[:,j] .- 1)) / Nz
    println(s)
    TV[j] = s
end
p = plot()
for j=1:length(t_list)
    if j<=5
        # dont plot t=0.216s
        plot!(p,C[:,j],zgrid, label="t=$(t_list[j])",dpi=300)
    end
end
#savefig(p, "./analy_1D_diffusion.png" )
# plot derivative approximation
# d1 = dkdz.(zgrid)
# d2 = zeros(Nz)
# for i=1:Nz
#     d2[i] = dkdz_approx(zgrid[i], 1 /(Nz-1))
# end
# p = plot(d1, zgrid, label="analy")
# plot!(d2, zgrid, label="approx")
# scatter!(d2, zgrid, label="")