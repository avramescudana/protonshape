using SpecialFunctions

"""
Gaussian radial envelope: exp(-alpha * r^2).
"""
function radial_envelope(r, alpha=0.05)
    return exp.(-alpha .* r.^2)
end

"""
Compute the angular modulation as a weighted sum of e^(im*phi) * J_n(kr).
"""
function angular_modulation(phi, r, coeff_dict, k=1.0)
    modulation = zeros(Float64, size(phi))
    for ((n, m), amp) in coeff_dict
        modulation .+= amp .* real(exp(im * m .* phi) .* besselj(n, k .* r))
    end
    return modulation
end

"""
Build the base density field as:
   density_base = envelope(r) * exp(angular_modulation(phi, r)).
"""
function build_base_density_2D_general(X, Y, coeff_dict; alpha=0.05, envelope_func=radial_envelope, k=1.0)
    r = sqrt.(X.^2 .+ Y.^2)
    phi = atan.(Y, X)
    phi = ifelse.(phi .< 0, phi .+ 2π, phi)

    env = envelope_func(r, alpha)
    ang_mod = angular_modulation(phi, r, coeff_dict, k)
    density_base = env .* exp.(ang_mod)
    return density_base
end
