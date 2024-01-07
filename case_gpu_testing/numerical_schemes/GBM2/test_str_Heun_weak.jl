μ=1.0; σ=1.0; S0=1.0; # specify the parameteres of the equivalent GBM to obtain analytical solution
Tend = 1
# This calculates the L1 error in the strong sense of the geometric Brownian motion
Nrep = 50000000
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
        S_star = S + 0.5*dt*S + S.*dW  # Stjnen(2002) equation(3.35)
        S = S + 0.5*dt*(0.5*S+0.5*S_star) + 0.5*(S+S_star).*dW 
    end   
    Smean = sum(S) / Nrep
    Sanaly = S0 * exp(μ*Tend)
    L1error = abs(Smean - Sanaly)
    println("N=$N, L1error=$(L1error)")
    L1errorlist[i] = L1error
end

p2 = plot(Nlist, L1errorlist,  xscale=:log10, yscale=:log10, minorgrid=true, label="",dpi=300)
scatter!(Nlist, L1errorlist, label="L1 error")
for i =1:length(Nlist) - 1
    tg = abs( (log10(L1errorlist[i+1]) - log10(L1errorlist[i]) ) / ( log10(Nlist[i+1]) - log10(Nlist[i]) )  )
    tg_text = @sprintf("%.3f", tg)
    annotate!(0.5*(Nlist[i]+Nlist[i+1]), 0.5*(L1errorlist[i]+L1errorlist[i+1]), tg_text)
end
plot!(p2, xlabel=L"N_T", ylabel=L"E_w(T)", title="L1 error at T=$(Tend)", titlefontsize=18, xguidefontsize=18, yguidefontsize=18)
savefig(p2, "./Heunweakerror.png" )
writedlm("./Heunrweakerror_result.txt", hcat(Nlist,L1errorlist),',')