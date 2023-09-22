using DifferentialEquations,Plots
α = 1
β = 1
u₀ = 1 / 2
f(u, p, t) = α * u
g(u, p, t) = β * u
dt = 1 // 2^(2)
tspan = (0.0, 1.0)
f_analytic(u₀, p, t, W) = u₀ * exp((α - (β^2) / 2) * t + β * W)
ff = SDEFunction(f, g, analytic = f_analytic)
prob = SDEProblem(ff, g, u₀, (0.0, 1.0))
sol = solve(prob, SRIW1(), adaptive=true)
plot(sol, plot_analytic=true)