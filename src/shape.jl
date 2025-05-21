using SpecialFunctions # For Bessel functions
using Roots  # For root-finding

"""
Gaussian radial envelope: exp(-α * r^2).
"""
function gaussenv(r, α=1.0)
    # return exp.(-α .* r.^2)
    return exp.(-r.^2 ./ (2*α)) ./ (2*π*α)
end
"""
Compute the n-th positive root of the Bessel function J_m.
"""
#TODO: optimize or tabulate the roots
function besseljzero(m, n; tol=1e-8)
    f(x) = besselj(m, x) 
    roots = []
    x_start = 0.1  
    step_size = 0.5 
    while length(roots) < n
        x_end = x_start + step_size  
        if f(x_start) * f(x_end) < 0  
            root = find_zero(f, (x_start, x_end), Bisection(), tol=tol)
            push!(roots, root)
            x_start = root + tol  
        else
            x_start = x_end  
        end
    end
    return roots[n]  
end

"""
Compute the weighted sum of e^(im*phi) * J_m(α_{mn}r/a).
"""
function circmemb(phi, r, coeff_dict, a=1.0)
    modulation = ones(size(phi))
    for ((m, n), amp) in coeff_dict
        α_mn = besseljzero(m, n)  
        modulation .+= amp .* real.(exp.(im * m .* phi) .* besselj.(m, α_mn ./ a .* r))
    end
    return modulation  
end

function circmemb_2D(X, Y, coeff_dict; a=1.0)
    r = sqrt.(X.^2 .+ Y.^2)
    phi = atan.(Y, X)
    phi = ifelse.(phi .< 0, phi .+ 2π, phi)
    memb = circmemb(phi, r, coeff_dict, a)
    return memb
end

"""
Build the thickness function Tp(phi, r) = gaussenv(r) * circmemb(phi, r).
"""
function Tp_2D(X, Y, coeff_dict; α=1.0, envfunc=gaussenv, a=1.0)
    r = sqrt.(X.^2 .+ Y.^2)
    phi = atan.(Y, X)
    phi = ifelse.(phi .< 0, phi .+ 2π, phi)

    env = envfunc(r, α)
    memb = circmemb(phi, r, coeff_dict, a)
    Tp = env .* memb # thickness function
    return Tp
end

"""
Shape thickness function used in the dipole cross section
"""

function Tp_shape(b, θb, p_shape)
    env = gaussenv(b, p_shape.α)
    memb = circmemb(θb, b, p_shape.coeff_dict, p_shape.a)
    Tp = env .* memb # thickness function
    return Tp
end