include("wavefunction.jl") # Overlap of photon and vector meson wavefunctions

function Agbw(r, b, z, Δ, p_wavefct, p_gbw)
	# Note that A(λ="L") for Q²=0

	# t1 = im * r * b
	# t1 has an overall factor of i
	t1 = r * b
	t2 = ΨᵥΨ(r, z, "T", p_wavefct)
	t3 = gbwdipole(r, b, p_gbw)
	t4 = besselj(0, b * Δ) * besselj(0, (1-z) * r * Δ)

	prod = t1 * t2 * t3 * t4
    return Symbolics.value(prod) # convert symbolic expression to numerical for MC integral
end

function diffractive(diff, dip, p_wavefct, p_gbw, p_mc)
    if diff == "coh"
        if dip == "GWB"
            include("src/gbwdipole.jl") # GWB dipole model

            xgbw = MCIntegration.Continuous(0,1)

            Δ_range = range(p_mc.Δmin, stop=p_mc.Δmax, length=p_mc.Δlen)
            t_range = Δ_range .* Δ_range 
            
            collect_int, collect_int_std = [], []
            
            for Δᵢ in Δ_range
                
                res = MCIntegration.integrate((xgbw, c)->Agbw(xgbw[1]*p_mc.rmax, xgbw[2]*p_mc.bmax, xgbw[3], Δᵢ, p_wavefct, p_gbw) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xgbw, dof = 3, solver=:vegas, neval=p_mc.neval)
            
                mean, std = res[1][1], res[1][2]

                push!(collect_int, mean)
                push!(collect_int_std, std)
            end
            
            units = (16/π) * 10^6 # [nb/GeV²]

            dσdt = abs.(collect_int .* collect_int) ./ units 
	        dσdt_err = 2 .* collect_int .* collect_int_std ./ units
        end
    end

    return t_range, dσdt, dσdt_err
end