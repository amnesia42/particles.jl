include("testeuler_lib.jl")
μ=0.5; σ=1.0; S0=1.0; # specify the parameteres of the GBM
Tend = 1
# This calculates the L1 error in the strong sense of the geometric Brownian motion
Nrep = 100000
Nlist = [2^i for i=3:9]
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
        S = S + dt*μ.*S  + σ*S.*dW + 0.5*σ^2*S .* (dW.*dW .- dt)
    end  
    Sanaly = S0 * ones(Nrep) .* exp.((μ - 0.5 * σ^2) * Tend .+ σ * W)
    L1error = sum(abs.(Sanaly- S )) / Nrep
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
plot!(p2, xlabel="N", ylabel=L"E|\hat{X}(T) - X(T)|", title="L1 error at T=$(Tend)", dpi=300)
savefig(p2, "./M1strongerror.png" )