include("testeuler_lib.jl")
μ=0.5; σ=1.0; S0=1.0; # specify the parameteres of the GBM
Tend = 1
# This calculates the L1 error in the strong sense of the geometric Brownian motion
Nrep = 5000000
p = plot()
Nlist = [2^i for i=3:6]
L1errorlist = zeros(length(Nlist))
for i in eachindex(Nlist)
    N = Nlist[i]
    dt = Tend/N
    S = zeros(Int(N)+1)
    Sanaly = zeros(Int(N)+1)
    for i=1:Nrep
        randnum_seq = randn(N)
        S .= S .+ get_gbm_euler(Tend, N, μ, σ, S0, randnum_seq)
        Sanaly .= Sanaly .+ get_gbm_analy(Tend, N, μ, σ, S0, randnum_seq)
    end
    S = S ./ Nrep
    Sanaly = Sanaly ./ Nrep
    #err = abs(S[end]-Sanaly[end]) #empirical analytical solution
    err = abs(S[end]-exp(μ*Tend)) # real analytical solution
    println("N=$N,err=$(err)")
    L1errorlist[i]=err
    plot!(p, range(0, Tend, Int(N)+1), S, label="euler, N=$(Int(N))")
end
plot!(p, range(0, Tend, 1000), exp.(μ .* range(0, Tend, 1000)), label="analy")

plot!(xlabel="t", ylabel=L"X_t")
savefig(p, "./multipleN.png" )

p2 = plot(Nlist, L1errorlist,  xscale=:log10, yscale=:log10, minorgrid=true, label="")
scatter!(Nlist, L1errorlist, label="L1 error")
for i =1:length(Nlist) - 1
    tg = abs( (log10(L1errorlist[i+1]) - log10(L1errorlist[i]) ) / ( log10(Nlist[i+1]) - log10(Nlist[i]) )  )
    tg_text = @sprintf("%.3f", tg)
    annotate!(0.5*(Nlist[i]+Nlist[i+1]), 0.5*(L1errorlist[i]+L1errorlist[i+1]), tg_text)
end
plot!(p2, xlabel="N", ylabel=L"|E[X_{T}] - E[X(T)]|", title="L1 error at t=$(Tend)")
savefig(p2, "./error.png" )
