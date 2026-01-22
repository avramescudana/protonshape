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
    # modulation = zeros(size(phi))
    for ((m, n), amp) in coeff_dict
        α_mn = besseljzero(m, n)  
        modulation .+= amp .* real.(exp.(im * m .* phi)) .* besselj.(m, α_mn ./ a .* r)
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

function Tp_shape(b, θb, p_shape; rotate=false)
    env = gaussenv(b, p_shape.α)
    if rotate
        θ_shift = 2π * rand()
        θb_rot = θb .+ θ_shift
        θb_rot = θb_rot .% (2π)
        memb = circmemb(θb_rot, b, p_shape.coeff_dict, p_shape.a)
    else
        memb = circmemb(θb, b, p_shape.coeff_dict, p_shape.a)
    end
    Tp = env .* memb
    return Tp
end

##############################################################################
# 3D SHAPE FLUCTUATIONS (SPHERICAL HARMONICS)
##############################################################################

"""
Compute angular modulation as a weighted sum of real spherical harmonics Y_l^m.
coeff_dict is a Dict with keys (l, m) and values being the amplitudes.
Uses SphericalHarmonics.jl for Y_l^m computation.
θ is the polar angle (colatitude, 0 to π), φ is the azimuthal angle (0 to 2π).
"""
function angular_modulation_3D(θ, φ, coeff_dict)
    modulation = zeros(size(θ))
    for ((l, m), amp) in coeff_dict
        if amp != 0
            # SphericalHarmonics.jl: sphericalharmonic(θ, φ, l, m) returns Y_l^m
            # Use RealHarmonics to get real spherical harmonics (matches Python's sph_harm.real)
            Ylm = real.(SphericalHarmonics.sphericalharmonic.(θ, φ, l, m, Ref(SphericalHarmonics.RealHarmonics())))
            modulation .+= amp .* Ylm
        end
    end
    return modulation
end

"""
3D angular membrane function using spherical harmonics (analog of circmemb for 2D).
Returns 1 + weighted sum of real spherical harmonics.
θ is polar angle (colatitude), φ is azimuthal angle.
Uses SphericalHarmonics.jl for Y_l^m computation.
"""
function sphermemb(θ, φ, coeff_dict)
    modulation = ones(size(θ))
    for ((l, m), amp) in coeff_dict
        if amp != 0
            Ylm = real.(SphericalHarmonics.sphericalharmonic.(θ, φ, l, m, Ref(SphericalHarmonics.RealHarmonics())))
            modulation .+= amp .* Ylm
        end
    end
    return modulation
end

"""
3D angular membrane function on Cartesian grid (analog of circmemb_2D).
X, Y, Z are 3D arrays of Cartesian coordinates.
"""
function sphermemb_3D(X, Y, Z, coeff_dict)
    r = sqrt.(X.^2 .+ Y.^2 .+ Z.^2)
    # Polar angle θ (colatitude): arccos(z/r)
    θ = similar(r)
    θ .= ifelse.(r .> 0, acos.(clamp.(Z ./ r, -1.0, 1.0)), 0.0)
    # Azimuthal angle φ
    φ = atan.(Y, X)
    φ = ifelse.(φ .< 0, φ .+ 2π, φ)

    memb = sphermemb(θ, φ, coeff_dict)
    return memb
end

"""
Generate a 3D noise field using a non-parametric continuous power spectrum.

Parameters:
  - shape: tuple (Nx, Ny, Nz) for grid shape
  - L: half-domain size (domain = [-L, L])
  - σ: scales the overall level of the power spectrum
  - seed: random seed (optional, uses global RNG if not provided)

Returns:
  - noise_field: real 3D noise array
"""
function generate_noise_field_3D(shape::Tuple{Int,Int,Int}, L::Real; σ::Real=2.0, seed::Union{Int,Nothing}=nothing)
    Nx, Ny, Nz = shape

    if seed !== nothing
        rng = Random.MersenneTwister(seed)
    else
        rng = Random.default_rng()
    end

    # FFT frequencies (wavenumbers) for each dimension
    # Grid spacing: dx = 2L/N, then k = 2π * fftfreq(N, dx)
    dx = 2*L / Nx
    dy = 2*L / Ny
    dz = 2*L / Nz
    kx = fftfreq(Nx, dx) .* 2π
    ky = fftfreq(Ny, dy) .* 2π
    kz = fftfreq(Nz, dz) .* 2π

    # Create 3D meshgrid of wavenumbers
    Kx = reshape(kx, Nx, 1, 1)
    Ky = reshape(ky, 1, Ny, 1)
    Kz = reshape(kz, 1, 1, Nz)

    Kmag = sqrt.(Kx.^2 .+ Ky.^2 .+ Kz.^2)
    k_max = maximum(Kmag)

    # Define a non-parametric continuous power spectrum with a peak
    k_peak = 4π
    # Build anchor points, filtering out any that exceed k_max
    k_candidates = [0.0, 0.5*k_peak, k_peak, 2*k_peak]
    S_candidates = [1*σ^2, 1*σ^2, 4.0*σ^2, 0.1*σ^2]

    # Keep only anchor points below k_max, then add k_max as final point
    mask = k_candidates .< k_max
    k_anchor = vcat(k_candidates[mask], k_max)
    S_anchor = vcat(S_candidates[mask], 1e-6*σ^2)

    # PCHIP interpolation in log space
    log_S_anchor = log.(S_anchor)
    itp = interpolate((k_anchor,), log_S_anchor, Gridded(Linear()))

    # Evaluate power spectrum at all k values
    S_continuous = exp.(itp.(Kmag)) .* exp.(-Kmag ./ k_peak)

    # Fourier filter: square root of S_continuous
    combined_filter = sqrt.(S_continuous)

    # Generate random complex noise
    noise_real = randn(rng, Nx, Ny, Nz)
    noise_imag = randn(rng, Nx, Ny, Nz)
    noise = noise_real .+ im .* noise_imag

    # Apply filter in Fourier space and transform back
    ρ_k = noise .* combined_filter
    noise_field = real.(ifft(ρ_k))

    return noise_field
end

"""
Helper function to compute FFT frequencies (like numpy.fft.fftfreq).
"""
function fftfreq(n::Int, d::Real=1.0)
    if iseven(n)
        return vcat(0:(n÷2-1), (-n÷2):-1) ./ (n * d)
    else
        return vcat(0:((n-1)÷2), (-(n-1)÷2):-1) ./ (n * d)
    end
end

"""
Combine base density with noise fluctuations:
    density = base_density * exp(amplitude * noise_field)
The exponent is clamped to avoid overflow.
"""
function combine_density_with_noise(base_density, noise_field; amplitude::Real=2.0, max_exponent::Real=20.0)
    exponent = amplitude .* noise_field
    exponent = clamp.(exponent, -max_exponent, max_exponent)
    return base_density .* exp.(exponent)
end

"""
Build 3D thickness function Tp(X, Y, Z) = gaussenv(r) * sphermemb(θ, φ).
This is the 3D analog of Tp_2D.

Parameters:
  - X, Y, Z: 3D arrays of Cartesian coordinates
  - coeff_dict: Dict with keys (l, m) and amplitude values
  - α: Gaussian width parameter
  - envfunc: radial envelope function (default: gaussenv)
  - add_noise: whether to add noise fluctuations
  - noise_σ: noise strength parameter
  - noise_amplitude: amplitude for combining noise
  - seed: random seed for noise
"""
function Tp_3D(X, Y, Z, coeff_dict; α::Real=1.0, envfunc=gaussenv,
               add_noise::Bool=false, noise_σ::Real=2.0, noise_amplitude::Real=2.0,
               seed::Union{Int,Nothing}=nothing)
    r = sqrt.(X.^2 .+ Y.^2 .+ Z.^2)
    # Polar angle θ (colatitude)
    θ = similar(r)
    θ .= ifelse.(r .> 0, acos.(clamp.(Z ./ r, -1.0, 1.0)), 0.0)
    # Azimuthal angle φ
    φ = atan.(Y, X)
    φ = ifelse.(φ .< 0, φ .+ 2π, φ)

    env = envfunc(r, α)
    memb = sphermemb(θ, φ, coeff_dict)
    Tp = env .* memb

    if add_noise
        shape = size(X)
        # Estimate L from grid (assume symmetric grid centered at origin)
        L = maximum(abs.(X))
        noise_field = generate_noise_field_3D(shape, L; σ=noise_σ, seed=seed)
        Tp = combine_density_with_noise(Tp, noise_field; amplitude=noise_amplitude)
    end

    return Tp
end

"""
3D shape thickness function used in the dipole cross section.
This is the 3D analog of Tp_shape.

Parameters:
  - r: radial distance
  - θ: polar angle (colatitude)
  - φ: azimuthal angle
  - p_shape: struct containing α and coeff_dict
  - rotate: if true, apply random rotation in φ
"""
function Tp_shape_3D(r, θ, φ, p_shape; rotate::Bool=false)
    env = gaussenv(r, p_shape.α)
    if rotate
        φ_shift = 2π * rand()
        φ_rot = mod.(φ .+ φ_shift, 2π)
        memb = sphermemb(θ, φ_rot, p_shape.coeff_dict)
    else
        memb = sphermemb(θ, φ, p_shape.coeff_dict)
    end
    Tp = env .* memb
    return Tp
end

"""
3D shape thickness function with Cartesian input.
Converts (x, y, z) to spherical coordinates and evaluates.

Parameters:
  - x, y, z: Cartesian coordinates (can be arrays)
  - p_shape: struct containing α and coeff_dict
  - rotate: if true, apply random rotation in φ
"""
function Tp_shape_3D_cartesian(x, y, z, p_shape; rotate::Bool=false)
    r = sqrt.(x.^2 .+ y.^2 .+ z.^2)
    θ = similar(r)
    θ .= ifelse.(r .> 0, acos.(clamp.(z ./ r, -1.0, 1.0)), 0.0)
    φ = atan.(y, x)
    φ = ifelse.(φ .< 0, φ .+ 2π, φ)

    return Tp_shape_3D(r, θ, φ, p_shape; rotate=rotate)
end

"""
Integrate the 3D density along the z-axis using the trapezoidal rule.
Returns a 2D array (projection onto the x-y plane).

Parameters:
  - density_3D: 3D array of density values (indexed as [x, y, z])
  - z_vals: 1D array of z coordinates
"""
function integrate_along_z(density_3D::Array{T,3}, z_vals::AbstractVector) where T
    Nx, Ny, Nz = size(density_3D)
    density_2D = zeros(T, Nx, Ny)

    # Trapezoidal rule integration along z
    for i in 1:Nx
        for j in 1:Ny
            density_2D[i, j] = trapz(z_vals, view(density_3D, i, j, :))
        end
    end

    return density_2D
end

"""
Simple trapezoidal integration.
"""
function trapz(x::AbstractVector, y::AbstractVector)
    n = length(x)
    if n < 2
        return zero(eltype(y))
    end
    s = zero(promote_type(eltype(x), eltype(y)))
    for i in 1:(n-1)
        s += 0.5 * (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return s
end