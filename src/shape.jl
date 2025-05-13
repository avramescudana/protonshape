using SpecialFunctions # For Bessel functions
using Roots  # For root-finding

"""
Gaussian radial envelope: exp(-alpha * r^2).
"""
function radial_envelope(r, alpha=0.05)
    return exp.(-alpha .* r.^2)
end

"""
Compute the n-th positive root of the Bessel function J_m.
"""
function besseljzero(m, n; tol=1e-8, max_iter=100)
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
function circular_membrane(phi, r, coeff_dict, a=1.0)
    modulation = ones(size(phi))
    for ((m, n), amp) in coeff_dict
        α_mn = besseljzero(m, n)  
        modulation .+= amp .* real.(exp.(im * m .* phi) .* besselj.(m, α_mn ./ a .* r))
    end
    return modulation  
end

"""
Build the base density field as:
   density_base = envelope(r) * circular_membrane(phi, r).
"""
function density_2D(X, Y, coeff_dict; alpha=0.05, envelope_func=radial_envelope, a=1.0)
    r = sqrt.(X.^2 .+ Y.^2)
    phi = atan.(Y, X)
    phi = ifelse.(phi .< 0, phi .+ 2π, phi)

    env = envelope_func(r, alpha)
    memb = circular_membrane(phi, r, coeff_dict, a)
    density_base = env .* memb
    return density_base
end
