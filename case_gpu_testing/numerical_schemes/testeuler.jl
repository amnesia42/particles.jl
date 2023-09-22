include("testeuler_lib.jl")
# This calculates the L1 error in the weak sense of the geometric Brownian motion

Tend = 1
Nrep = 50000
Nlist = [2^i for i=5:10]
#Nlist = [10,100,1000,10000,100000]
L1errorlist = zeros(length(Nlist))
μ=0.5; σ=1.0; S0=1.0; # specify the parameteres of the GBM


for i in eachindex(Nlist)
    L1errorlist[i] = get_gbmmeanerror_L1norm(Tend, Nlist[i],μ, σ, S0, Nrep)
end

println(L1errorlist)

p2 = plot(Nlist, L1errorlist,  xscale=:log10, yscale=:log10, minorgrid=true, label="")
scatter!(Nlist, L1errorlist, label="L1 error")
#xlims!(0.5e+1, 1.5e+6)
#ylims!(1e-6, 5e-2)
plot!(p2, xlabel="N", ylabel=L"|E(X_{t_N}) - E(X_T)|", title="L1 error at t=$(Tend)")
# calculate the tangent of the curve
# add the annotated absolute tangent value
for i =1:length(Nlist) - 1
    tg = abs( (log10(L1errorlist[i+1]) - log10(L1errorlist[i]) ) / ( log10(Nlist[i+1]) - log10(Nlist[i]) )  )
    tg_text = @sprintf("%.3f", tg)
    annotate!(0.3*(Nlist[i]+Nlist[i+1]), 0.3*(L1errorlist[i]+L1errorlist[i+1]), tg_text)
end

savefig(p2, "./L1meanerror_vs_N_Nrep=$(Nrep).png" )
