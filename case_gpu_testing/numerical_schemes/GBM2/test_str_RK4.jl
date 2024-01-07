μ=1.0; σ=1.0; S0=1.0; # specify the parameteres of the equivalent GBM to obtain analytical solution
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
        K0 = 0.5*S
        G0 = S
        Y0 = S + 0.5*dt*K0 + 0.5*G0.*dW
        K1 = 0.5*Y0
        G1 = Y0
        Y1 = S + 0.5*dt*K1 + 0.5*G1.*dW
        K2 = 0.5*Y1
        G2 = Y1
        Y2 = S + dt*K2 + G2.*dW
        K3 = 0.5*Y2
        G3 = Y2
        
        S = S + dt/6*(K0+2*K1+2*K2+K3) + 1.0/6*(G0+2*G1+2*G2+G3).*dW
    end   
    Sanaly = S0 * ones(Nrep) .* exp.((μ - 0.5 * σ^2) * Tend .+ σ * W)
    L1error = sum(abs.(Sanaly- S )) / Nrep
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
plot!(p2, xlabel=L"N_T", ylabel=L"E_s(T)", title="L1 error at T=$(Tend)", titlefontsize=18, xguidefontsize=18, yguidefontsize=18)
savefig(p2, "./RK4strongerror.png" )