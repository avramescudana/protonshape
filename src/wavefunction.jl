"""
Functions
"""

function ϕ(r, z, λ, p)
	if λ=="T"
		norm = p.Nₜ
	elseif λ=="L"
		norm = p.Nₗ
	end

	factor = norm * z * (1-z)
	term1 = - p.mf * p.mf * p.R² / (8 * z * (1-z))
	term2 = - 2 * z * (1-z) * r * r / p.R²
	term3 = p.mf * p.mf * p.R² / 2
	term_exp = exp(term1 + term2 + term3)
    
	return factor * term_exp
end

ϵ(z, p) = sqrt(z * (1-z) * p.Q² + p.mf * p.mf)

# Symbolic derivatives
∂ᵣ = Differential(r)
∂ᵣϕₜ(r, z, p) = expand_derivatives(∂ᵣ(ϕ(r, z, "T", p)))
∇ᵣϕₗ(r, z, p) = 1/r * expand_derivatives(∂ᵣ(ϕ(r, z, "L", p))) + expand_derivatives(∂ᵣ(∂ᵣ(ϕ(r, z, "L", p))))

function ΨᵥΨ(r, z, λ, p)
	factor = p.êf * p.e * p.Nc / π
    
	if λ=="T"
		factorλ = 1 / (z * (1-z))
		term1 = p.mf * p.mf * besselk(0, ϵ(z, p) * r) * ϕ(r, z, λ, p)
		term2 = - (z * z + (1-z) * (1-z)) * ϵ(z, p) * besselk(1, ϵ(z, p) * r) * ∂ᵣϕₜ(r, z, p)
	elseif λ=="L"
		factorλ = 2 * sqrt(p.Q²) * z * (1-z) * besselk(0, ϵ(z, p) * r)
		term1 = p.Mᵥ * ϕ(r, z, λ, p)
		term2 = p.δ / (p.Mᵥ * z * (1-z)) * (p.mf * p.mf * ϕ(r, z, λ, p) - ∇ᵣϕₗ(r, z, p))
	end

	return factor * factorλ * (term1 + term2) 
end