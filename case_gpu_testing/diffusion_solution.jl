include("ConcentrationCalculationLibrary.jl")
using LegendrePolynomials
using Plots
z0=0.5
Nz = 101
#t_list = [0.036 0.072 0.108 0.1728]
t_list = [0.036]
zgrid = range(0,1,Nz)

# plot analytical solution 
C = ones(Nz, length(t_list))
for j = 1:length(t_list)
    t = t_list[j]
    for n=1:15
        for i=1:Nz
            δC= (2*n+1)*Pl(2*zgrid[i]-1, n)*Pl(2*z0-1,n)*exp(-6*n*(n+1)*t)
            C[i,j] += δC
            if i==101
                println("n={$n}, C(0.5)=$(δC)")
            end
        end
    end
end
p = plot(C,zgrid)

# plot derivative approximation
# d1 = dkdz.(zgrid)
# d2 = zeros(Nz)
# for i=1:Nz
#     d2[i] = dkdz_approx(zgrid[i], 1 /(Nz-1))
# end
# p = plot(d1, zgrid, label="analy")
# plot!(d2, zgrid, label="approx")
# scatter!(d2, zgrid, label="")