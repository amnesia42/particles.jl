using CUDA, Test, BenchmarkTools
using LinearAlgebra, SpecialFunctions
using DelimitedFiles
using Plots
using PyCall, Printf
np = pyimport("numpy")

# diffusion functions not in a hurry

# kernel estimator using the reflection method to have order-1 consistency
function Kernel_Estimator(x, Locparticles, bandwidth)
    # kernel methods
    function Gauss_kernel(u; d=1)
        return (2*pi)^(-0.5*d) * exp(-0.5*u*u)
    end

    function Epa_kernel(u; d=1)
        nu_d = 2*pi^(0.5*d) / d / gamma(0.5*d) 
        if u*u > 1
            return 0.
        else
            return 0.5 / nu_d * (d+2) * (1-u*u)
        end
    end

    #Kernel = Gauss_kernel
    Kernel = Epa_kernel

    d = 1
    Nparticles = length(Locparticles)
    p = 0
    if x<bandwidth
        #print("reflection applies.")
        for i in eachindex(Locparticles)
            u = (Locparticles[i] - x)/bandwidth
            p += Kernel(u)
            # the kernel span outside the boundary; max when the sample is at the boundary
            # aka the deflection length
            d_reflection = bandwidth-(Locparticles[i]-0) 
            if x<d_reflection
                u_reflection = (Locparticles[i] + x)/bandwidth
                p+= Kernel(u_reflection)
            end
        end
        temp = p/(Nparticles*bandwidth^d)
    elseif 1-x< bandwidth
        #print("reflection applies.")
        for i in eachindex(Locparticles)
            u = (Locparticles[i] - x)/bandwidth
            p += Kernel(u)
            d_reflection = bandwidth-(1 -Locparticles[i])
            if 1-x<d_reflection
                u_reflection = (2-x - Locparticles[i])/bandwidth
                p+= Kernel(u_reflection)
            end
        end
        temp = p/(Nparticles*bandwidth^d)
    else
        for i in eachindex(Locparticles)
            u = (Locparticles[i] - x)/bandwidth
            p+= Kernel(u)
        end
        temp = p/(Nparticles*bandwidth^d)
    end
    return temp
end