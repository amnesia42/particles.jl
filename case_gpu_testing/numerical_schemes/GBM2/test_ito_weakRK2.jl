μ=1.0; σ=1.0; S0=1.0; # specify the parameteres of the GBM
function f(x)
    return μ * x
end
function g(x)
    return σ * x
end
Tend = 1
# This calculates the L1 error in the strong sense of the geometric Brownian motion
Nrep = 1000
Nlist = [2^i for i=3:6]
L1errorlist = zeros(length(Nlist))
for i in eachindex(Nlist)
    N = Nlist[i]
    dt = Tend/N
    S = S0 .* ones(Nrep)
    W = zeros(Nrep)
    t = 0
    for j=1:N
        dW = randn(rng, Nrep) * sqrt(dt)
        W = W + dW
        S_tilde = S + f.(S) * dt + g.(S) .* dW # Saerkka(2019)) equation(8.98)
        S_plus = S + f.(S) * dt + g.(S) * sqrt(dt)
        S_minus = S + f.(S) * dt - g.(S) * sqrt(dt)
        S = S + 0.5*(f.(S)+f.(S_tilde))*dt + 0.25*(g.(S_plus)+2*g.(S)+g.(S_minus)).*dW
              + 0.25/sqrt(dt)*(g.(S_plus)-g.(S_minus)).*(dW.*dW .- dt)
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