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
    prod = t1 * t2 * t3 * t4 * t5
	if part == "real"
		return Symbolics.value(real(prod))
	elseif part == "imag"
		return Symbolics.value(imag(prod))
	end
end

function threaded_loop(run_threads, range, body)
    if run_threads
        @threads for i in range
            body(i)
        end
    else
        for i in range
            body(i)
        end
    end
end

function diffractive(diff, dip, p_wavefct, p_gbw, p_cq, p_mc; run_threads=false)
    if dip == "GWB"
        #TODO: rewrite this part to use the same structure as CQ
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

        if diff == "coh"
            dσdt = abs.(collect_int .* collect_int) ./ units 
            dσdt_mcerr = 2 .* collect_int .* collect_int_std ./ units
        elseif diff == "incoh"
            println("Incoherent cross section not implemented for GWB dipole model")
        end
    elseif dip == "CQ"
        xqc = MCIntegration.Continuous(0,1)

        Δ_range = range(p_mc.Δmin, stop=p_mc.Δmax, length=p_mc.Δlen)
        t_range = Δ_range .* Δ_range 

        collect_abs2 = [Float64[] for _ in 1:p_cq.Nsamples]
        collect_A = [ComplexF64[] for _ in 1:length(Δ_range)]

        threaded_loop(run_threads, 1:p_cq.Nsamples, ibq -> begin
            bqc_sample = sample_bqc(p_cq)  
            abs2_for_sample = Float64[]  

            for (i, Δᵢ) in enumerate(Δ_range)
                resqc_re = MCIntegration.integrate((xqc, c) -> Aqc(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, bqc_sample, xqc[4], Δᵢ, Tp, p_wavefct, p_gbw, p_cq, "real") * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

                resqc_imag = MCIntegration.integrate((xqc, c) -> Aqc(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, bqc_sample, xqc[4], Δᵢ, Tp, p_wavefct, p_gbw, p_cq, "imag") * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

                A_sample = (-resqc_imag[1][1] + resqc_re[1][1] * im) / 2.0  # A ∝ i e^(-iB)

                push!(abs2_for_sample, abs2(A_sample))  
                push!(collect_A[i], A_sample)
            end

            collect_abs2[ibq] = abs2_for_sample
        end)

        abs2_mean = [abs2(mean(A_samples)) for A_samples in collect_A]
        mean_abs2 = [mean(abs2_for_t) for abs2_for_t in eachcol(collect_abs2)][1]
        factor = 389.38 / (16π) # GeV^-2 = 389.38 mb, overall 1(16π) factor

        std_A = [std(abs2.(A_samples)) for A_samples in collect_A]
        N = length(collect_A[1])  
        error_abs2_mean = [std / sqrt(N) for std in std_A]

        dσdt_coh = abs2_mean .* factor
        dσdt_coh_err = error_abs2_mean .* factor
        dσdt_incoh = (mean_abs2 -  abs2_mean) .* factor
    end
    if diff == "coh" || diff == "incoh"
        # TODO: implement output only for coh or incoh
        return t_range
    elseif diff == "coh+incoh"
        return t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh
    end
end