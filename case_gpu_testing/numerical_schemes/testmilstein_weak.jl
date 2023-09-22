include("testeuler_lib.jl")
μ=0.5; σ=1.0; S0=1.0; # specify the parameteres of the GBM
Tend = 1
# This calculates the L1 error in the strong sense of the geometric Brownian motion
Nrep = 50000000
Nlist = [2^i for i=3:6]
L1errorlist = zeros(length(Nlist))
for i in eachindex(Nlist)
    N = Nlist[i]
    dt = Tend/N
    S = S0 .* ones(Nrep)
    for j=1:N
        dW = randn(rng, Nrep) * sqrt(dt)
        S = S + dt*μ.*S  + σ*S.*dW + 0.5*σ^2*S .* (dW.*dW .- dt)
    end   
    L1error = abs(exp(μ*Tend) - sum(S)/Nrep  )
    println("N=$N, L1error=$(L1error), S_mean=$(sum(S)/Nrep)")
    L1errorlist[i] = L1error
end

p2 = plot(Nlist, L1errorlist,  xscale=:log10, yscale=:log10, minorgrid=true, label="", dpi=300)
scatter!(Nlist, L1errorlist, label="L1 error")
for i =1:length(Nlist) - 1
    tg = abs( (log10(L1errorlist[i+1]) - log10(L1errorlist[i]) ) / ( log10(Nlist[i+1]) - log10(Nlist[i]) )  )
    tg_text = @sprintf("%.3f", tg)
    annotate!(0.5*(Nlist[i]+Nlist[i+1]), 0.5*(L1errorlist[i]+L1errorlist[i+1]), tg_text)
end
plot!(p2, xlabel="N", ylabel=L"|E[\hat{X}(T)] - E[X(T)]|", title="L1 error at T=$(Tend)")
savefig(p2, "./M1weakerror.png")