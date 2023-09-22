μ=1.0; σ=1.0; S0=1.0; # specify the parameteres of the GBM
Tend = 1
# This calculates the L1 error in the strong sense of the geometric Brownian motion
Nrep = 1000
Nlist = [2^i for i=4:9]
L1errorlist = zeros(length(Nlist))
for i in eachindex(Nlist)
    N = Nlist[i]
    dt = Tend/N
    S = S0 .* ones(Nrep)
    W = zeros(Nrep)
    t = 0
    for j=1:N
        dW = randn(rng, Nrep) * sqrt(dt)
        dZ = randn(rng, Nrep) * sqrt(dt) # another sequence of random number required for RK2
        W = W + dW
        GAMMA = 1.0/sqrt(2)*S.*dZ + dt/4*S + 0.5*S.*dZ.*dZ # Graewe(2011)) equation(21-22) Heemink(2018) equation(9.58)
        S = S + sqrt(2)*dW.*GAMMA + 1.0/sqrt(2)*S.*(dZ-dW) + 0.5*GAMMA*dt + (0.5*S+1.0*GAMMA).*dZ.*dZ + S.*(dZ-0.5*dW).^2
    end   
    Smean = sum(S) / length(S)
    Sanaly = exp(μ*Tend)
    L1error = sum(abs.(Sanaly- Smean)) / Nrep
    println("N=$N, L1error=$(L1error)")
    L1errorlist[i] = L1error
end

p2 = plot(Nlist, L1errorlist,  xscale=:log10, yscale=:log10, minorgrid=true, label="")
scatter!(Nlist, L1errorlist, label="L1 error")
for i =1:length(Nlist) - 1
    tg = abs( (log10(L1errorlist[i+1]) - log10(L1errorlist[i]) ) / ( log10(Nlist[i+1]) - log10(Nlist[i]) )  )
    tg_text = @sprintf("%.3f", tg)
    annotate!(0.5*(Nlist[i]+Nlist[i+1]), 0.5*(L1errorlist[i]+L1errorlist[i+1]), tg_text)
end
plot!(p2, xlabel="N", ylabel=L"|E[X_{T}] - E[X(T)]|", title="L1 error at t=$(Tend)")
savefig(p2, "./Rk2weakerror.svg" )