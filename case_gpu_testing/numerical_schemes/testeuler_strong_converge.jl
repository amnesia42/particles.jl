include("testeuler_lib.jl")

Tend = 1
NLength = 100
# This calculates the L1 error in the strong sense of the geometric Brownian motion
Nrep = 100000
Nlist = [2^i for i=3:9]
L1errorlist = zeros(length(Nlist))
for i in eachindex(Nlist)
    N = Nlist[i]
    L1error = 0
    for i=1:Nrep
        dW_list = randn(rng, N)
        L1error += get_gbmerror_L1norm(1, N, 0.5, 1.0, 1.0, dW_list)
    end
    L1error = L1error / Nrep
    println("N=$N, L1error=$(L1error)")
    L1errorlist[i] = L1error
end

p2 = plot(Nlist, L1errorlist,  xscale=:log10, yscale=:log10, minorgrid=true, label="", dpi=300)
scatter!(Nlist, L1errorlist, label="L1 error")
plot!(p2, xlabel="N", ylabel=L"E[|X(T)-̂\hat{X}(T)]", title="L1 error at T=$(Tend)")
# calculate the tangent of the curve
# add the annotated absolute tangent value
for i =1:length(Nlist) - 1
    tg = abs( (log10(L1errorlist[i+1]) - log10(L1errorlist[i]) ) / ( log10(Nlist[i+1]) - log10(Nlist[i]) )  )
    tg_text = @sprintf("%.3f", tg)
    annotate!(0.5*(Nlist[i]+Nlist[i+1]), 0.5*(L1errorlist[i]+L1errorlist[i+1]), tg_text)
end
savefig(p2, "./eulerstrongerror.png" )
