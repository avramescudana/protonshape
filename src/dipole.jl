# include("shape.jl") # Circular membrane model 
using Random

function Qₛ(p)
	return (p.x₀ / p.xₚ) ^ (p.Λ/2) # [GeV²]
end

"""
GBW dipole model
"""

function T(b, p)
	norm = 1 / (2π * p.Bₚ)
	term_exp = exp(- b * b / (2 * p.Bₚ))
	return norm * term_exp
end

function gbwdipole(r, b, p)
	term_exp = 1 - exp(- p.N₀ * r * r * Qₛ(p) * Qₛ(p) * T(b, p) / 4)
	# return p.σ₀ * term_exp
    return 2 * term_exp
end

"""
Constituent quark model
"""

function sample_bqc(p)
	# (bx, by) for each quark i according to a Gaussian
	σqc = √p.Bqc
	gauss_qc = Normal(0, σqc)
	bqc = [(rand(gauss_qc), rand(gauss_qc)) for _ in 1:p.Nq]
	return bqc
end

function Tq(b, p)
    norm_squared = b[1]^2 + b[2]^2
    return (1 / (2π * p.Bq)) * exp(-norm_squared / (2 * p.Bq))
end

function Tp(b, θb, bqc, p)
	bx = b * cos(θb)
    by = b * sin(θb)
    vecb = (bx, by)
	
    sum = 0.0
    for bᵢ in bqc
        delta = (vecb[1] - bᵢ[1], vecb[2] - bᵢ[2])
        sum += Tq(delta, p)
    end
    return sum / p.Nq
end

function compute_Tp_grid(bx_vals, by_vals, bqc, p)

    Tp_vals = Matrix{Float64}(undef, length(by_vals), length(bx_vals))

    for (i, by) in enumerate(by_vals), (j, bx) in enumerate(bx_vals)
        b = hypot(bx, by)
        θ = atan(by, bx)
        Tp_vals[i, j] = Tp(b, θ, bqc, p)
    end

    return Tp_vals
end

function cqdipole(r, b, θb, bqc, Tp, p_cq)
	term_Tp = Tp(b, θb, bqc, p_cq)
    term_exp = 1 - exp(- p_cq.N₀ * r * r * term_Tp)
    return 2 * term_exp
end

"""
Shape thickness function
"""

function shapedipole(r, b, θb, Tp_shape, p_shape)
	term_Tp = Tp_shape(b, θb, p_shape; p_shape.rotate)
    term_exp = 1 - exp(- p_shape.N₀ * r * r * term_Tp)
    return 2 * term_exp
end

function sample_amp_dict_same_mn(p)
    damp = Normal(0, p.σ)
	amps = rand(damp, p.Nsamples)
    # coeff_dict = copy(p.coeff_dict)
    # coeff_dict[p.mn] = amps
    dicts = [Dict(p.mn => amp) for amp in amps]
	# return coeff_dict
    return dicts
end

function sample_amp_dict_samem_multin(p)
    damp = Normal(0, p.σ)
    dicts = []
    for nsample in 1:p.Nsamples
        amp_dict = Dict()
        for n in 1:p.nvals
            amp_dict[(p.mn[1], n)] = rand(damp)
        end
        push!(dicts, amp_dict)
    end
    return dicts
end