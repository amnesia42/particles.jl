include("testeuler_lib.jl")

Tend = 1
N = 1000
μ=0.5; σ=1.0; S0=1.0; # specify the parameteres of the GBM
grid_ratio = 10

randnum = randn(rng, N)
S_analy = get_gbm_analy(Tend, N, μ, σ, S0, randnum)
S_euler = get_gbm_euler(Tend, N, μ, σ, S0, randnum; fine_to_coarse_grid_ratio=grid_ratio)

p = plot(range(0, Tend, N+1), S_analy, label="analy")
plot!(p,range(0, Tend, Int(N/grid_ratio)+1), S_euler, label="euler")
plot!(p, xlabel="t", ylabel=L"X_t", size=(400,300), dpi=300)
savefig(p, "euler_versus_analy_adjustsize.png")


