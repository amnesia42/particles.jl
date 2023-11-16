using CUDA, Adapt, Test, BenchmarkTools
using LinearAlgebra, SpecialFunctions
using DelimitedFiles
using Plots
using Printf
using StatsBase

# diffusion functions not in a hurry
function k(z)
    # parabola diffusion
    return 6*z.*(1 .- z)
end

function dkdz(z)
    return 6*(1 .- 2*z)
end

function dkdz_approx(z,h)
    if z+h>1
        return (1.5*k(z)-2*k(z-h)+0.5*k(z-2h))/h
    elseif z-h<0
        return (-1.5*k(z)+2*k(z+h)-0.5*k(z+2*h))/h
    else
        return (k(z+h)-k(z-h))/(2*h)
    end 
end

function get_next_location(p, dt; h=0.002)
    # diffusion functions not in a hurry
    function k(z)
        # parabola diffusion
        if 0≤z≤1
            return 6*z*(1 - z)
        else
            return 0.
        end
    end

    function dkdz_approx(z,h)
        if z+h>1
            return (1.5*k(z)-2*k(z-h)+0.5*k(z-2h))/h
        elseif z-h<0
            return (-1.5*k(z)+2*k(z+h)-0.5*k(z+2*h))/h
        elseif z>1 || z<0
            return 0.
        else
            return (k(z+h)-k(z-h))/(2*h)
        end 
    end

    # Numerical schemes
    # dX = f*dt + σ*dW
    function get_Δx_euler(p)
        # Graewe 2011, eq.9
        return dkdz_approx(p,h)*dt + sqrt(2*k(p))*(randn()*sqrt(dt))
    end
    function get_Δx_m1(p)
        # Graewe 2011, eq. 11, Spivaskovskaya 2007 eq.34
        dW = randn()*sqrt(dt)
        f = dkdz_approx(p,h)
        σ = sqrt(2*k(p))
        σtimesσprime = dkdz_approx(p,h)
        return f*dt + σ*dW + 0.5*σtimesσprime*(dW*dW-dt)
    end
    function get_Δx_heun(p)
        #Stjnen PhD Thesis 2002, eq.3.35
        dW = randn()*sqrt(dt)  # \tilde f = (f-0.5*σ*σ') = 0.5*dkdz
        σ = sqrt(2)*sqrt(k(p))    #scheme="heun"
        ftilde = 0.5*dkdz_approx(p, h)
        pnext = p + ftilde * dt + σ*dW 
        return 0.5*(ftilde+0.5*dkdz_approx(pnext,h))*dt + 0.5*(σ + sqrt(2*k(pnext)) ) *dW
    end
    # function get_Δx_heun(p)
    #     dW = randn()*sqrt(dt)  # \tilde f = (f-0.5*σ*σ') = 0.5*dkdz
    #     f = dkdz_approx(p,h)
    #     σ = sqrt(2*k(p))    #scheme="heun"
    #     σprime = dkdz_approx(p,h)/σ
    #     fnext = f - 0.5*σ*σprime
    #     pnext = p + fnext * dt + σ*dW
        
    #     fnext = dkdz_approx(pnext,h)
    #     σnext = sqrt(2*k(pnext))
    #     σprimenext = dkdz_approx(pnext,h)/σnext
    #     fnextnext = fnext - 0.5*σnext*σprimenext
    #     return fnextnext*dt + σnext*dW
    # end
    function get_Δx_RK4(p)
        # Stijnen PhD Thesis 2002, eq. 3.37
        # under Stratonovich representation
        # σtimesσprime = dkdz_approx(p, h)
        # ftilde = f-0.5*σtimesσprime = 0.5*dkdz_approx(p, h); 
        dW = randn()*sqrt(dt)
        K0 = 0.5*dkdz_approx(p,h)
        G0 = sqrt(2*k(p))
        pnext0 = p + 0.5*K0*dt + 0.5*G0*dW
        K1 = 0.5*dkdz_approx(pnext0,h)
        G1 = sqrt(2*k(pnext0))
        pnext1 = p + 0.5*K1*dt + 0.5*G1*dW
        K2 = 0.5*dkdz_approx(pnext1,h)
        G2 = sqrt(2*k(pnext1))
        pnext2 = p + K2*dt + G2*dW
        K3 = 0.5*dkdz_approx(pnext2,h)
        G3 = sqrt(2*k(pnext2))
        return 1.0/6*(K0+2*K1+2*K2+K3)*dt +1.0/6*(G0+2*G1+2*G2+G3)*dW
    end

    #if scheme=="euler"
    dp = get_Δx_euler(p)
    #elseif scheme=="m1"
    #dp = get_Δx_m1(p)
    #elseif scheme == "heun"
    #dp = get_Δx_heun(p)
    #elseif scheme == "RK4"
    #dp = get_Δx_RK4(p)
    #end
    if p+dp>1 || p+dp<0
        p = get_next_location(p, 0.5*dt)
        p = get_next_location(p, 0.5*dt)
    else
        return p+dp
    end
end

function get_bandwidth(Locparticles)
    Nparticles = length(Locparticles)
    x_avg = sum(Locparticles)/length(Locparticles)
    x_delta = Locparticles .- x_avg
    std_sample = sqrt(1/(Nparticles - 1)* dot(x_delta, x_delta))
    bandwidth = std_sample * Nparticles^(-0.2)
    return bandwidth
end

function Box_Estimator(xgrid, Locparticles)
    hist = fit(Histogram, Locparticles, xgrid)
    return hist.weights ./ length(Locparticles)
end

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

# function new_Kernel_Estimator(xgrid, cgrid, Locparticle, bandwidth)
#     # kernel methods
#     function Gauss_kernel(u; d=1)
#         return (2*pi)^(-0.5*d) * exp(-0.5*u*u)
#     end

#     function Epa_kernel(u; d=1)
#         nu_d = 2*pi^(0.5*d) / d / gamma(0.5*d) 
#         if u*u > 1
#             return 0.
#         else
#             return 0.5 / nu_d * (d+2) * (1-u*u)
#         end
#     end

#     #Kernel = Gauss_kernel
#     Kernel = Epa_kernel
#     #cgrid = CUDA.zeros(length(xgrid))
#     d = 1
#     for i in eachindex(xgrid)
#         x = xgrid[i]
#         if x<bandwidth
#             #print("reflection applies.")
#             u = (Locparticle - x)/bandwidth
#             cgrid[i] += Kernel(u)
#             # the kernel span outside the boundary; max when the sample is at the boundary
#             # aka the deflection length
#             d_reflection = bandwidth - (Locparticle-0) 
#             if x<d_reflection
#                 u_reflection = (Locparticle + x)/bandwidth
#                 cgrid[i]+= Kernel(u_reflection)
#             end
#         elseif 1-x< bandwidth
#             #print("reflection applies.")
#             u = (Locparticle - x)/bandwidth
#             cgrid[i] += Kernel(u)
#             d_reflection = bandwidth-(1 -Locparticle)
#             if 1-x<d_reflection
#                 u_reflection = (2-x - Locparticle)/bandwidth
#                 cgrid[i]+= Kernel(u_reflection)
#             end
#         else
#             u = (Locparticle - x)/bandwidth
#             cgrid[i]+= Kernel(u)
#         end
#     end
#     return cgrid
# end

function new_Kernel_Estimator(x, Locparticle, bandwidth)
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
    #cgrid = CUDA.zeros(length(xgrid))
    d = 1
    p = 0 

    if x<bandwidth
        #print("reflection applies.")
        u = (Locparticle - x)/bandwidth
        p += Kernel(u)
        # the kernel span outside the boundary; max when the sample is at the boundary
        # aka the deflection length
        d_reflection = bandwidth - (Locparticle-0) 
        if x<d_reflection
            u_reflection = (Locparticle + x)/bandwidth
            p += Kernel(u_reflection)
        end
    elseif 1-x< bandwidth
        #print("reflection applies.")
        u = (Locparticle - x)/bandwidth
        p += Kernel(u)
        d_reflection = bandwidth-(1 -Locparticle)
        if 1-x<d_reflection
            u_reflection = (2-x - Locparticle)/bandwidth
            p += Kernel(u_reflection)
        end
    else
        u = (Locparticle - x)/bandwidth
        p += Kernel(u)
    end
    return p
end