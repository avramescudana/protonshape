include("wavefunction.jl") # Overlap of photon and vector meson wavefunctions
include("dipole.jl") # GWB and CQ dipole models

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

function Aqc(r, b, θb, bqc, z, Δ, Tp, p_wavefct, p_gbw, p_cq, part)
	# Note that A(λ="L") for Q²=0
	
	# t1 = im * r * b
	# t1 has an overall factor of i/2
	t1 = r * b
	t2 = ΨᵥΨ(r, z, "T", p_wavefct)
	t3 = cqdipole(r, b, θb, bqc, Tp, p_gbw, p_cq)
	t4 = exp(- im * b * Δ * cos(θb))
	t5 = besselj(0, (1-z) * r * Δ)
	if part == "real"
		return Symbolics.value(real(t1 * t2 * t3 * t4 * t5))
	elseif part == "imag"
		return Symbolics.value(imag(t1 * t2 * t3 * t4 * t5))
	end
end

function diffractive(diff, dip, p_wavefct, p_gbw, p_cq, p_mc)
    if diff == "coh"
        if dip == "GWB"
            xgbw = MCIntegration.Continuous(0,1)

            Δ_range = range(p_mc.Δmin, stop=p_mc.Δmax, length=p_mc.Δlen)
            t_range = Δ_range .* Δ_range 
            
            collect_int, collect_int_std = [], []
            
            for Δᵢ in Δ_range
                
                res = MCIntegration.integrate((xgbw, c)->Agbw(xgbw[1]*p_mc.rmax, xgbw[2]*p_mc.bmax, xgbw[3], Δᵢ, p_wavefct, p_gbw) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xgbw, dof = 3, solver=:vegas, neval=p_mc.neval, niters=p_mc.niters)
            
                mean, std = res[1][1], res[1][2]

                push!(collect_int, mean)
                push!(collect_int_std, std)
            end
            
            units = (16/π) * 10^6 # [nb/GeV²]

            dσdt = abs.(collect_int .* collect_int) ./ units 
	        dσdt_err = 2 .* collect_int .* collect_int_std ./ units
        elseif dip == "CQ"
            xqc = MCIntegration.Continuous(0,1)
	
            Δ_range = range(p_mc.Δmin, stop=p_mc.Δmax, length=p_mc.Δlen)
            t_range = Δ_range .* Δ_range 
            
            collect_int, collect_int_std = [], []
            
            for Δᵢ in Δ_range
                bqc_samples = [sample_bqc(p_cq) for _ in 1:p_cq.Nsamples]
                int_samples = ComplexF64[]
            
                for bqc_sample in bqc_samples
                
                    resqc_re = MCIntegration.integrate((xqc, c)->Aqc(xqc[1]*p_mc.rmax, xqc[2]*p_mc.bmax, xqc[3]*p_mc.θbmax, bqc_sample, xqc[4], Δᵢ, Tp, p_wavefct, p_gbw, p_cq, "real") * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver=:vegas, neval=p_mc.neval, niters=p_mc.niters) 

                    resqc_imag = MCIntegration.integrate((xqc, c)->Aqc(xqc[1]*p_mc.rmax, xqc[2]*p_mc.bmax, xqc[3]*p_mc.θbmax, bqc_sample, xqc[4], Δᵢ, Tp, p_wavefct, p_gbw, p_cq, "imag") * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver=:vegas, neval=p_mc.neval, niters=p_mc.niters)
                
                    meanqc = resqc_re[1][1] + resqc_imag[1][1] * im
            
                    push!(int_samples, meanqc)
                end

                mean_abs2 = mean(abs2, int_samples)
                abs2_mean = abs2(mean(int_samples))
            
                var_qc = mean_abs2 - abs2_mean
                std_qc = sqrt(var_qc)

                # push!(collect_var_qc, var_qc)
                # push!(collect_std_qc, std_qc)

                push!(collect_int, abs2_mean)
                push!(collect_int_std, std_qc)

            units = (16/π) * 10^6 # [nb/GeV²]

            dσdt = abs.(collect_int .* collect_int) ./ units 
            dσdt_err = 2 .* collect_int .* collect_int_std ./ units
            end
        end
    end

    return t_range, dσdt, dσdt_err
end