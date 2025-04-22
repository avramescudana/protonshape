"""
Functions
"""

function Qₛ(p)
	return (p.x₀ / p.xₚ) ^ (p.Λ/2) # [GeV²]
end

function T(b, p)
	norm = 1 / (2π * p.Bₚ)
	term_exp = exp(- b * b / (2 * p.Bₚ))
	return norm * term_exp
end

function gbwdipole(r, b, p)
	term_exp = 1 - exp(- p.N₀ * r * r * Qₛ(p) * Qₛ(p) * T(b, p) / 4)
	return p.σ₀ * term_exp
end

