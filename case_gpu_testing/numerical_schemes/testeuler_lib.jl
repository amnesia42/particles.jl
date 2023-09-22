using DelimitedFiles
using Random
using Plots, LaTeXStrings, Printf

rng = MersenneTwister(123)


function get_gbm_analy(T, N, μ, σ, S0, randnum_seq; fine_to_coarse_grid_ratio=1)
    if N == length(randnum_seq)
        dt = T / N 
        Ncoarse = Int(N/fine_to_coarse_grid_ratio)
        S = zeros(Ncoarse+1)
        S[1] = S0
        t = 0.0
        W = 0.0
        for i in 1:Ncoarse
            t += dt*fine_to_coarse_grid_ratio
            start_index = Int(1 + (i-1) * fine_to_coarse_grid_ratio)
            end_index = Int(start_index + fine_to_coarse_grid_ratio - 1)
            dW = sum(randnum_seq[start_index:end_index]) * sqrt(dt) # Brownian motion increment
            W += dW
            S[i+1] = S0 * exp((μ - 0.5 * σ^2) * t + σ * W)
        end
        return S
    else
        throw("The length of the given random number sequence doesn't match the input N. ")
    end
end

function get_gbm_euler(T, N, μ, σ, S0, randnum_seq; fine_to_coarse_grid_ratio = 1)
    if N == length(randnum_seq)
        dt = T / N      # the dt that matches the size of the randnum_seq
        NLength_euler = Int(N/fine_to_coarse_grid_ratio)
        S = zeros(NLength_euler+1)
        S[1] = S0
        for i = 1:NLength_euler
            start_index = 1 + (i-1) * fine_to_coarse_grid_ratio
            end_index = start_index + fine_to_coarse_grid_ratio - 1
            dW = sum(randnum_seq[start_index:end_index]) * sqrt(dt) # Brownian motion increment
            S[i+1] = S[i] + μ*S[i]*(dt*fine_to_coarse_grid_ratio) +σ*S[i]*dW
        end
        return S
    else
        throw("The length of the given random number sequence doesn't match the input N. ")
    end
end

function get_gbmerror_L1norm(Tend, NLength, μ, σ, S0, randnum_seq)
    xgrid = zeros(NLength+1)
    xnumgrid = zeros(NLength+1)
    xgrid .+= get_gbm_analy(Tend, NLength, μ, σ, S0, randnum_seq)
    xnumgrid .+= get_gbm_euler(Tend, NLength, μ, σ, S0, randnum_seq)
    # compute L1 norm at t=Tend
    L1_norm = abs(xgrid[end] - xnumgrid[end])
    return L1_norm
end

function get_gbmmeanerror_L1norm(Tend, NLength, μ, σ, S0, Nrep)
    dt = Tend / NLength          
    S_sum = zeros(NLength+1)
    for j in 1:Nrep
        randnum_seq = randn(NLength)
        S_sum = S_sum .+ get_gbm_euler(Tend,NLength,μ,σ,S0,randnum_seq)
    end
    S_sum = S_sum ./ Nrep
    println(S_sum[end])
    return abs(S_sum[end] - S0*exp(μ*Tend))
end