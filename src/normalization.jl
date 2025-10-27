using Optim
using Interpolations

# Compute χ² for a given N₀
function optimal_norm_N0(N₀_values, params_shape, params_wavefct, params_mc, params_run, coherent_data_path)
    best_χsq = Inf
    best_N₀ = nothing
    χsq_list = Float64[]
    N₀_list = Float64[]

    isdir(params_run.savepath) || mkpath(params_run.savepath)
    outdir_path = params_run.savepath * "/" * params_run.outdir
    isdir(outdir_path) || mkpath(outdir_path)

    for N₀_val in N₀_values
        params_shape_norm = merge(params_shape, (N₀ = N₀_val,))
        Δ_val = 0.0
        params_mc_Δ_val = merge(params_mc, (Δmin = 0.0, Δmax = 0.0, Δlen = 1))

        diffractive("coh+incoh", "shapeamp", params_wavefct, params_mc_Δ_val; p_shape=params_shape_norm, p_run=params_run, Δ_values=[Δ_val])

        Nsamples = params_shape.Nsamples

        t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err = compute_cross_sections(outdir_path, [Δ_val], Nsamples, params_run, Nsamples)

        tcent_coh_hera, dσdt_coh_hera, Δtot_coh_hera = read_coherent_data(coherent_data_path)

        model_val = dσdt_coh[1]
        data_val = dσdt_coh_hera[1]
        data_err = Δtot_coh_hera[1]

        χsq = ((model_val - data_val) / data_err)^2

        push!(χsq_list, χsq)
        push!(N₀_list, N₀_val)

        if χsq < best_χsq
            best_χsq = χsq
            best_N₀ = N₀_val
        end
    end

    return best_N₀, best_χsq, N₀_list, χsq_list
end

# Compute χ² at Δ=0 for a given N₀
function chisq_for_N₀_at_Δ₀(N₀, params_shape, params_wavefct, params_mc, params_run, params_norm; Δ₀=0.0)
    try
        params_shape_norm = merge(params_shape, (N₀ = N₀,))
        params_mc_Δ0 = merge(params_mc, (Δmin = Δ₀, Δmax = Δ₀, Δlen = 1))

        p_run_iter = if params_norm.unique_outdirs
            outdir = joinpath(params_run.outdir, "N0_$(round(N₀; digits=6))")
            merge(params_run, (outdir = outdir,))
        else
            params_run
        end

        diffractive("coh+incoh", "shapeamp", params_wavefct, params_mc_Δ0; p_shape=params_shape_norm, p_run=p_run_iter)

        compute_dir = joinpath(p_run_iter.savepath, p_run_iter.outdir)
        endswith(compute_dir, "/") && (compute_dir = compute_dir[1:end-1])

        Nsamples = params_shape.Nsamples
        t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err = compute_cross_sections(compute_dir, [Δ₀], Nsamples, p_run_iter, Nsamples)

        tcent_coh_hera, dσdt_coh_hera, Δtot_coh_hera = read_coherent_data(params_norm.coherent_data_path)

        model_val = dσdt_coh[1]
        data_val  = dσdt_coh_hera[1]
        data_err  = Δtot_coh_hera[1]

        χsq = ((model_val - data_val) / data_err)^2

        return isfinite(χsq) ? χsq : Inf
    catch
        return Inf
    end
end

# Adaptive search starting at N₀=1.0 that can go lower or higher
function find_best_N₀_at_Δ₀_adaptive(params_shape, params_wavefct, params_mc, params_run, params_norm)

    cache = Dict{Float64,Float64}()
    obj(N) = get!(cache, N) do
        chisq_for_N₀_at_Δ₀(N, params_shape, params_wavefct, params_mc, params_run, params_norm)
    end

    x₀ = clamp(params_norm.start, params_norm.min_N₀, params_norm.max_N₀)
    f₀ = obj(x₀)

    xᵣ, fᵣ = x₀, f₀
    n = 0
    while n < params_norm.max_expansions
        xtry = clamp(xᵣ * params_norm.step_factor, params_norm.min_N₀, params_norm.max_N₀)
        if xtry == xᵣ; break; end
        ftry = obj(xtry)
        if ftry < fᵣ
            xᵣ, fᵣ = xtry, ftry
            n += 1
        else
            break
        end
    end

    xₗ, fₗ = x₀, f₀
    n = 0
    while n < params_norm.max_expansions
        xtry = clamp(xₗ / params_norm.step_factor, params_norm.min_N₀, params_norm.max_N₀)
        if xtry == xₗ; break; end
        ftry = obj(xtry)
        if ftry < fₗ
            xₗ, fₗ = xtry, ftry
            n += 1
        else
            break
        end
    end

    a = min(xₗ, xᵣ); b = max(xₗ, xᵣ)
    if a == b
        a = max(params_norm.min_N₀, x₀/params_norm.step_factor)
        b = min(params_norm.max_N₀, x₀*params_norm.step_factor)
    end
    println("Optimizing in interval: a = $a, b = $b")
    if a >= b
        error("x_lower must be less than x_upper: a = $a, b = $b")
    end

    result = optimize(obj, a, b, Brent(); rel_tol=params_norm.brent_reltol)
    best_N₀ = Optim.minimizer(result)
    best_χsq = Optim.minimum(result)

    # return best_N₀, best_χsq, (a,b), cache
     return best_N₀
end