"""
Main functions
"""

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

function A(r, b, θb, z, Δ, Tp, p_wavefct, dip, part; 
            bqc=nothing, p_cq=nothing, p_shape=nothing)
    t1 = r * b
    t2 = ΨᵥΨ(r, z, "T", p_wavefct)
    if dip == "CQ"
        t3 = cqdipole(r, b, θb, bqc, Tp, p_cq)
    elseif dip == "shape"
        t3 = shapedipole(r, b, θb, Tp, p_shape)
    else
        error("Unknown dipole type: $dip")
    end
    t4 = exp(-im * b * Δ * cos(θb))
    t5 = besselj(0, (1-z) * r * Δ)
    prod = t1 * t2 * t3 * t4 * t5
    if part == "real"
        return Symbolics.value(real(prod))
    elseif part == "imag"
        return Symbolics.value(imag(prod))
    else
        error("Unknown part: $part")
    end
end

function diffractive(diff, dip, p_wavefct, p_mc; p_gbw=nothing, p_cq=nothing, p_shape=nothing, p_run=nothing, Δ_values=nothing)
    if dip == "GWB"
        #TODO: rewrite this part to use the same structure as CQ
        xgbw = MCIntegration.Continuous(0,1)

        if Δ_values === nothing
            Δ_range = range(p_mc.Δmin, stop=p_mc.Δmax, length=p_mc.Δlen)
        else
            Δ_range = Δ_values
        end
        t_range = Δ_range .* Δ_range 
        
        collect_int, collect_int_std = [], []
        
        #TODO: use threads for this part
        @showprogress for Δᵢ in Δ_range
            res = suppress_mcoutput(() -> begin
                MCIntegration.integrate((xgbw, c)->Agbw(xgbw[1]*p_mc.rmax, xgbw[2]*p_mc.bmax, xgbw[3], Δᵢ, p_wavefct, p_gbw) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xgbw, dof = 3, solver=:vegas, neval=p_mc.neval, niters=p_mc.niters)
            end, redirect_output = true, output_file = "mcintegrate_log.txt")


            mean, std = res[1][1], res[1][2]
            mean, std = mean * π, std * π

            push!(collect_int, mean)
            push!(collect_int_std, std)
        end
        
        factor = 389.38 / (16π) # GeV^-2 = 389.38 mb, overall 1(16π) factor 

        if diff == "coh"
            dσdt = abs.(collect_int .* collect_int) .* factor 
            dσdt_mcerr = 2 .* collect_int .* collect_int_std .* factor
        elseif diff == "incoh"
            println("Incoherent cross section not implemented for GWB dipole model")
        end
        return t_range, dσdt, dσdt_mcerr
    elseif dip == "CQ" || dip == "shape" || dip == "shapeamp"
        xqc = MCIntegration.Continuous(0,1)

        Δ_range = range(p_mc.Δmin, stop=p_mc.Δmax, length=p_mc.Δlen)
        t_range = Δ_range .* Δ_range 

        if dip == "CQ"
            collect_abs2 = [Float64[] for _ in 1:p_cq.Nsamples]
            collect_A = [ComplexF64[] for _ in 1:length(Δ_range)]

            threaded_loop(p_run.run_threads, 1:p_cq.Nsamples, ibq -> begin
                bqc_sample = sample_bqc(p_cq)  
                abs2_for_sample = Float64[]  

                for (i, Δᵢ) in enumerate(Δ_range)
                    # log_file_re = "mcintegrate_re_log_thread_$(ibq)_index_$(i).log"
                    # log_file_im = "mcintegrate_im_log_thread_$(ibq)_index_$(i).log"
                    
                    # resqc_re = suppress_mcoutput_thread(() -> begin
                    #     MCIntegration.integrate((xqc, c) -> Aqc(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, bqc_sample, xqc[4], Δᵢ, Tp, p_wavefct, p_gbw, p_cq, "real") * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)
                    # end, redirect_output = true)

                    # resqc_imag = suppress_mcoutput_thread(() -> begin
                    #     MCIntegration.integrate((xqc, c) -> Aqc(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, bqc_sample, xqc[4], Δᵢ, Tp, p_wavefct, p_gbw, p_cq, "imag") * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)
                    # end, redirect_output = true)

                    resqc_re = MCIntegration.integrate((xqc, c) -> A(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, xqc[4], Δᵢ, Tp, p_wavefct, "CQ", "real"; bqc=bqc_sample, p_cq=p_cq) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

                    resqc_imag = MCIntegration.integrate((xqc, c) -> A(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, xqc[4], Δᵢ, Tp, p_wavefct, "CQ", "imag"; bqc=bqc_sample, p_cq=p_cq) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

                    A_sample = (-resqc_imag[1][1] + resqc_re[1][1] * im) / 2.0  # A ∝ i e^(-iB)

                    push!(abs2_for_sample, abs2(A_sample))  
                    push!(collect_A[i], A_sample)
                end

                collect_abs2[ibq] = abs2_for_sample
            end)
        elseif dip == "shape"
            collect_abs2 = Vector{Float64}(undef, length(Δ_range))
            collect_A = [ComplexF64[] for _ in 1:length(Δ_range)]

            threaded_loop(p_run.run_threads, 1:length(Δ_range), i -> begin
                Δᵢ = Δ_range[i]

                resqc_re = MCIntegration.integrate((xqc, c) -> A(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, xqc[4], Δᵢ, Tp_shape, p_wavefct, "shape", "real"; p_shape=p_shape) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

                resqc_imag = MCIntegration.integrate((xqc, c) -> A(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, xqc[4], Δᵢ, Tp_shape, p_wavefct, "shape", "imag"; p_shape=p_shape) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

                A_sample = (-resqc_imag[1][1] + resqc_re[1][1] * im) / 2.0 # A ∝ i e^(-iB)

                collect_abs2[i] = abs2(A_sample)
                push!(collect_A[i], A_sample)
            end)
        elseif dip == "shapeamp"
            if p_run.run == "local"
                collect_abs2 = [Float64[] for _ in 1:p_shape.Nsamples]
                collect_A = [ComplexF64[] for _ in 1:length(Δ_range)]

                if p_shape.type=="samemn"
                    coeff_dicts = sample_amp_dict_same_mn(p_shape)
                elseif p_shape.type=="samem_multin"
                    coeff_dicts = sample_amp_dict_samem_multin(p_shape)
                else
                    error("Unknown sampling type: $(p_shape.type)")
                end

                threaded_loop(p_run.run_threads, 1:p_shape.Nsamples, iamp -> begin 
                    if p_run.run == "local"
                        abs2_for_sample = Float64[]  
                    end

                    coeff_dict_amp = merge(p_shape, (coeff_dict=coeff_dicts[iamp],))

                    for (i, Δᵢ) in enumerate(Δ_range)
        
                        resqc_re = MCIntegration.integrate((xqc, c) -> A(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, xqc[4], Δᵢ, Tp_shape, p_wavefct, "shape", "real"; p_shape=coeff_dict_amp) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

                        resqc_imag = MCIntegration.integrate((xqc, c) -> A(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, xqc[4], Δᵢ, Tp_shape, p_wavefct, "shape", "imag"; p_shape=coeff_dict_amp) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

                        A_sample = (-resqc_imag[1][1] + resqc_re[1][1] * im) / 2.0  # A ∝ i e^(-iB)

                        push!(abs2_for_sample, abs2(A_sample))  
                        push!(collect_A[i], A_sample)

                        if p_run.savefile
                           isdir(p_run.savepath) || mkpath(p_run.savepath)
                            outdir_path = p_run.savepath * "/" * p_run.outdir
                            isdir(outdir_path) || mkpath(outdir_path)
                            filename = joinpath(outdir_path, "A_delta$(i)_config$(iamp).jld2")
                            @info "Saved A_sample" filename=filename Δ=Δᵢ iamp=iamp
                            @save filename A_sample Δᵢ iamp
                        end
                    end

                    collect_abs2[iamp] = abs2_for_sample
                end)
            elseif p_run.run == "remote"
                coeff_dicts = p_run.amp_dict
                iamp = p_run.arrayindex
                # coeff_dict_amp = merge(p_shape, (coeff_dict=coeff_dicts[iamp],))
                # Accept either:
                # - coeff_dicts::AbstractVector where coeff_dicts[i] is a Dict for config i
                # - coeff_dicts::Dict{Tuple{Int,Int},Float64} representing a single config used for all indices
                if coeff_dicts isa AbstractVector
                    coeff_dict_amp = merge(p_shape, (coeff_dict=coeff_dicts[iamp],))
                elseif coeff_dicts isa Dict
                    # single dictionary of (m,n)=>amp provided; use it for this job
                    coeff_dict_amp = merge(p_shape, (coeff_dict=coeff_dicts,))
                else
                    error("Unexpected type for p_run.amp_dict: $(typeof(coeff_dicts))")
                end

                for (i, Δᵢ) in enumerate(Δ_range)
    
                    resqc_re = MCIntegration.integrate((xqc, c) -> A(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, xqc[4], Δᵢ, Tp_shape, p_wavefct, "shape", "real"; p_shape=coeff_dict_amp) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

                    resqc_imag = MCIntegration.integrate((xqc, c) -> A(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, xqc[4], Δᵢ, Tp_shape, p_wavefct, "shape", "imag"; p_shape=coeff_dict_amp) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

                    A_sample = (-resqc_imag[1][1] + resqc_re[1][1] * im) / 2.0  # A ∝ i e^(-iB)

                    isdir(p_run.savepath) || mkpath(p_run.savepath)
                    outdir_path = p_run.savepath * "/" * p_run.outdir
                    isdir(outdir_path) || mkpath(outdir_path)
                    filename = joinpath(outdir_path, "A_delta$(i)_config$(iamp).jld2")
                    @info "Saved A_sample" filename=filename Δ=Δᵢ iamp=iamp
                    # if p_run.jobtype=="single"
                        # filename = joinpath(outdir_path, "A_delta$(i)_sample$(iamp).jld2")
                    # elseif p_run.jobtype=="array"
                        # filename = joinpath(outdir_path, "A_delta$(i)_sample$(p_run.arrayindex).jld2")
                    # else
                    #     error("Unknown job type: $(p_run.jobtype)")
                    # end
                    @save filename A_sample Δᵢ iamp
                    # params_file = joinpath(outdir_path, "params_config$(iamp).jld2")
                    # @save params_file diff dip p_wavefct p_mc p_gbw p_cq p_shape p_run
                    # open(params_file, "w") do io
                    #     serialize(io, (diff, dip, p_wavefct, p_mc, p_gbw, p_cq, p_shape, p_run))
                    # end
                end
            end
        end

        if p_run.run == "local"
            abs2_mean = [abs2(mean(A_samples)) for A_samples in collect_A]
            factor = 389.38 / (16π) # GeV^-2 = 389.38 mb, overall 1(16π) factor

            if diff == "coh" 
                mean_A = [mean(A_samples) for A_samples in collect_A]
                sem_A = [std(A_samples) / sqrt(length(A_samples)) for A_samples in collect_A]

                dσdt_coh = abs2.(mean_A) .* factor
                dσdt_coh_err = abs.(2 .* mean_A .* sem_A .* factor)

                return t_range, dσdt_coh, dσdt_coh_err

            elseif diff == "incoh"
                mean_abs2 = [mean(abs2_for_t) for abs2_for_t in eachcol(collect_abs2)][1]
                dσdt_incoh = (mean_abs2 -  abs2_mean) .* factor

                collect_abs2_matrix = hcat(collect_abs2...)
                mean_abs2_err = [std(collect_abs2_matrix[i, :]) / sqrt(length(collect_abs2_matrix[i, :])) for i in eachindex(collect_abs2_matrix[:, 1])]
                dσdt_incoh_err = sqrt.((mean_abs2_err .* factor) .^2 .+ (error_abs2_mean .* factor) .^2)

                return t_range, dσdt_incoh, dσdt_incoh_err

            elseif diff == "coh+incoh"
                mean_A = [mean(A_samples) for A_samples in collect_A]
                sem_A = [std(A_samples) / sqrt(length(A_samples)) for A_samples in collect_A]

                dσdt_coh = abs2.(mean_A) .* factor
                dσdt_coh_err = abs.(2 .* mean_A .* sem_A .* factor)
                # dσdt_coh_err = [std(abs2.(A_samples)) / sqrt(length(A_samples)) * factor for A_samples in collect_A]

                mean_abs2 = [mean(abs2_for_t) for abs2_for_t in eachcol(collect_abs2)][1]
                dσdt_incoh = (mean_abs2 -  abs2_mean) .* factor

                # std_A = [std(abs2.(A_samples)) for A_samples in collect_A]
                # error_abs2_mean = [std / sqrt(length(collect_A[1])  ) for std in std_A]
                # collect_abs2_matrix = hcat(collect_abs2...)
                # mean_abs2_err = [std(collect_abs2_matrix[i, :]) / sqrt(length(collect_abs2_matrix[i, :])) for i in eachindex(collect_abs2_matrix[:, 1])]
                # dσdt_incoh_err = sqrt.((mean_abs2_err .* factor) .^2 .+ (error_abs2_mean .* factor) .^2)

                dσdt_incoh_err = compute_incoh_error(collect_abs2, collect_A, factor; error_method=p_mc.error_method)

                # Half-sample variance 
                # collect_abs2_matrix = hcat(collect_abs2...) 

                # nconf = size(collect_abs2_matrix, 2)
                # half = div(nconf, 2)

                # var_full  = [var(collect_abs2_matrix[i, :]) for i in eachindex(collect_abs2_matrix[:, 1])]
                # var_half1 = [var(collect_abs2_matrix[i, 1:half]) for i in eachindex(collect_abs2_matrix[:, 1])]
                # var_half2 = [var(collect_abs2_matrix[i, half+1:end]) for i in eachindex(collect_abs2_matrix[:, 1])]

                # err1 = abs.(var_half1 .- var_full)
                # err2 = abs.(var_half2 .- var_full)
                # dσdt_incoh_err = 0.5 .* (err1 .+ err2) .* factor

                return t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err 
            end
        end
    end



    # elseif p_run.run == "remote"
    #     #TODO: add all dip modes
    #     if dip == "shapeamp"
    #         collect_A = [ComplexF64[] for _ in 1:length(Δ_range)]

    #         if p_shape.type=="samemn"
    #             coeff_dicts = sample_amp_dict_same_mn(p_shape)
    #         elseif p_shape.type=="samem_multin"
    #             coeff_dicts = sample_amp_dict_samem_multin(p_shape)
    #         else
    #             error("Unknown sampling type: $(p_shape.type)")
    #         end

    #         for iamp in range(p_shape.Nsamples)
    #             coeff_dict_amp = merge(p_shape, (coeff_dict=coeff_dicts[iamp],))

    #             for (i, Δᵢ) in enumerate(Δ_range)
    
    #                 resqc_re = MCIntegration.integrate((xqc, c) -> A(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, xqc[4], Δᵢ, Tp_shape, p_wavefct, "shape", "real"; p_shape=coeff_dict_amp) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

    #                 resqc_imag = MCIntegration.integrate((xqc, c) -> A(xqc[1] * p_mc.rmax, xqc[2] * p_mc.bmax, xqc[3] * p_mc.θbmax, xqc[4], Δᵢ, Tp_shape, p_wavefct, "shape", "imag"; p_shape=coeff_dict_amp) * (p_mc.rmax^2) * (p_mc.bmax^2); var = xqc, dof = 4, solver = :vegas, neval = p_mc.neval, niters = p_mc.niters)

    #                 A_sample = (-resqc_imag[1][1] + resqc_re[1][1] * im) / 2.0  # A ∝ i e^(-iB)

    #                 if p_run.savefile
    #                     outdir_path = p_run.savepath + p_run.outdir
    #                     isdir(outdir_path) || mkpath(outdir_path)
    #                     filename = joinpath(outdir, "A_delta$(lpad(i,4,'0'))_sample$(lpad(iamp,4,'0')).jld2")
    #                     @save filename A_sample Δᵢ iamp
    #                 else
    #                     push!(collect_A[i], A_sample)
    #                 end
    #             end

    #             if p_run.savefile == false
    #                 return t_range, collect_A
    #             end
    #         end
    #     end
    # end
    
end

"""
Helper functions
"""

function threaded_loop(run_threads, range, body)
    if run_threads
        # progress = Progress(length(range), desc="Progress")

        # @showprogress @threads for i in range
        @threads for i in range
            body(i)
            # lock(output_lock) do
            #     next!(progress)  
            # end
        end

        # finish!(progress)
    else
        @showprogress for i in range
            body(i)
        end
    end
end

function suppress_mcoutput(body; redirect_output = true, output_file = nothing)
    if redirect_output
        file = output_file === nothing ? "/dev/null" : output_file

        open(file, "w") do io
            redirect_stdout(io) do
                redirect_stderr(io) do
                    return body()
                end
            end
        end
    else
        return body()
    end
end

function suppress_mcoutput_thread(body; redirect_output = true)
    if redirect_output
        # temp_file = tempname()
        temp_file = "/dev/null"

        try
            open(temp_file, "w") do io
                redirect_stdout(io) do
                    redirect_stderr(io) do
                        return body()
                    end
                end
            end
        finally
            # isfile(temp_file) && rm(temp_file)
        end
    else
        return body()
    end
end

function compute_incoh_error(collect_abs2, collect_A, factor; error_method="standard")
    # collect_abs2: list of arrays, each array is abs2(A) for a config
    # collect_A: list of arrays, each array is A for a config

    collect_abs2_matrix = hcat(collect_abs2...)
    nconf = size(collect_abs2_matrix, 2)
    half = div(nconf, 2)

    if error_method == "standard"
        # Standard error of the mean
        std_A = [std(abs2.(A_samples)) for A_samples in collect_A]
        error_abs2_mean = [std / sqrt(length(collect_A[1])  ) for std in std_A]
        mean_abs2_err = [std(collect_abs2_matrix[i, :]) / sqrt(length(collect_abs2_matrix[i, :])) for i in eachindex(collect_abs2_matrix[:, 1])]
        dσdt_incoh_err = sqrt.((mean_abs2_err .* factor) .^2 .+ (error_abs2_mean .* factor) .^2)
        return dσdt_incoh_err

    elseif error_method == "halfsample"
        # Half-sample variance
        var_full  = [var(collect_abs2_matrix[i, :]) for i in eachindex(collect_abs2_matrix[:, 1])]
        var_half1 = [var(collect_abs2_matrix[i, 1:half]) for i in eachindex(collect_abs2_matrix[:, 1])]
        var_half2 = [var(collect_abs2_matrix[i, half+1:end]) for i in eachindex(collect_abs2_matrix[:, 1])]
        err1 = abs.(var_half1 .- var_full)
        err2 = abs.(var_half2 .- var_full)
        dσdt_incoh_err = 0.5 .* (err1 .+ err2) .* factor
        return dσdt_incoh_err

    elseif error_method == "jackknife"
        # Jackknife error estimation
        n = nconf
        jk_means = [mean(deleteat!(collect_abs2_matrix[i, :], j)) for i in eachindex(collect_abs2_matrix[:, 1]), j in 1:n]
        jk_mean = [mean(jk_means[i, :]) for i in eachindex(collect_abs2_matrix[:, 1])]
        jk_var = [(n-1)/n * sum((jk_means[i, :] .- jk_mean[i]).^2) for i in eachindex(collect_abs2_matrix[:, 1])]
        dσdt_incoh_err = sqrt.(jk_var) .* factor
        return dσdt_incoh_err

    else
        error("Unknown error_method: $error_method")
    end
end

# function compute_cross_sections(outdir::String, Δ_range::AbstractVector, Nsamples::Int, p_run::NamedTuple, nconfigs::Int, multiple_configs::Int)
#     NΔ = length(Δ_range)
#     collect_A = [ComplexF64[] for _ in 1:NΔ]

#     config_dirs = multiple_configs == 1 ? filter(f -> isdir(joinpath(outdir, f)), readdir(outdir)) : [""]
    
#     for config in config_dirs
#         config_path = multiple_configs == 1 ? joinpath(outdir, config) : outdir

#         for i in 1:NΔ
#             for iamp in 1:nconfigs
#                 filename = joinpath(config_path, "A_delta$(i)_sample$(iamp).jld2")

#                 if isfile(filename)
#                     data = JLD2.load(filename)
#                     A_sample = data["A_sample"]
#                     push!(collect_A[i], A_sample)
#                 else
#                     @warn "Missing file: $filename"
#                 end
#             end
#         end
#     end

#     factor = 389.38 / (16π)
#     t_range = Δ_range .^ 2

#     mean_A = [mean(A_samples) for A_samples in collect_A]
#     sem_A = [std(A_samples) / sqrt(length(A_samples)) for A_samples in collect_A]
#     abs2_mean = abs2.(mean_A)

#     dσdt_coh = abs2_mean .* factor
#     dσdt_coh_err = abs.(2 .* mean_A .* sem_A .* factor)

#     mean_abs2 = [mean(abs2.(A_samples)) for A_samples in collect_A]
#     dσdt_incoh = (mean_abs2 .- abs2_mean) .* factor
#     dσdt_incoh_err = [std(abs2.(A_samples)) / sqrt(length(A_samples)) * factor for A_samples in collect_A]

#     return t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err
# end

function compute_cross_sections(outdir::String, Δ_range::AbstractVector, Nsamples::Int, p_run::NamedTuple, nconfigs::Int)
    NΔ = length(Δ_range)
    collect_A = [ComplexF64[] for _ in 1:NΔ]

    for i in 1:NΔ
        for iamp in 1:nconfigs
            filename = joinpath(outdir, "A_delta$(i)_config$(iamp).jld2")

            if isfile(filename)
                data = JLD2.load(filename)
                A_sample = data["A_sample"]
                push!(collect_A[i], A_sample)
            else
                @warn "Missing file: $filename"
            end
        end
    end

    factor = 389.38 / (16π)
    t_range = Δ_range .^ 2

    mean_A = [mean(A_samples) for A_samples in collect_A]
    sem_A = [std(A_samples) / sqrt(length(A_samples)) for A_samples in collect_A]
    abs2_mean = abs2.(mean_A)

    dσdt_coh = abs2_mean .* factor
    dσdt_coh_err = abs.(2 .* mean_A .* sem_A .* factor)

    mean_abs2 = [mean(abs2.(A_samples)) for A_samples in collect_A]
    dσdt_incoh = (mean_abs2 .- abs2_mean) .* factor
    dσdt_incoh_err = [std(abs2.(A_samples)) / sqrt(length(A_samples)) * factor for A_samples in collect_A]

    println("Computed cross sections for path: $outdir")

    return t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err
end