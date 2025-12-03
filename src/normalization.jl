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
function χsq_for_N₀_at_Δ₀(N₀, params_shape, params_wavefct, params_mc, params_run, params_norm)
    try
        Δ₀ = params_norm.Δ₀
        params_shape_norm = merge(params_shape, (N₀ = N₀, Nsamples = params_norm.nsamples_norm))
        params_mc_Δ₀ = merge(params_mc, (Δmin = Δ₀, Δmax = Δ₀, Δlen = 1))

        p_run_iter = if params_norm.unique_outdirs
            outdir = joinpath(params_run.outdir, "N0_$(round(N₀; digits=6))")
            merge(params_run, (outdir = outdir,))
        else
            params_run
        end

        p_run_iter_local = merge(p_run_iter, (run = "local", savefile = true,
    run_threads = false,))

        diffractive("coh+incoh", "shapeamp", params_wavefct, params_mc_Δ₀; p_shape=params_shape_norm, p_run=p_run_iter_local)

        compute_dir = joinpath(p_run_iter.savepath, p_run_iter.outdir)
        endswith(compute_dir, "/") && (compute_dir = compute_dir[1:end-1])

        Nsamples = params_norm.nsamples_norm
        t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err = compute_cross_sections(compute_dir, [Δ₀], Nsamples, p_run_iter_local, Nsamples)

        tcent_coh_hera, dσdt_coh_hera, Δtot_coh_hera = read_coherent_data(params_norm.coherent_data_path)

        model_val = dσdt_coh[1]
        data_val  = dσdt_coh_hera[1]
        data_err  = Δtot_coh_hera[1]

        χsq = ((model_val - data_val) / data_err)^2

        return isfinite(χsq) ? χsq : Inf
    catch e
        @warn "χsq_for_N₀_at_Δ₀ failed" N₀=N₀ exception=e
        return Inf
    end
end

# Adaptive search starting at N₀=1.0 that can go lower or higher
function find_best_N₀_at_Δ₀_adaptive(params_shape, params_wavefct, params_mc, params_run, params_norm)

    cache = Dict{Float64,Float64}()

    # Obj prints every evaluation so you can trace what the optimizer is asking for.
    obj(N) = get!(cache, N) do
        println("[eval] χ² at N₀ = $(N)")
        f = χsq_for_N₀_at_Δ₀(N, params_shape, params_wavefct, params_mc, params_run, params_norm)
        println("[eval] result: χ²($(N)) = $(f)")
        f
    end

    x₀ = clamp(params_norm.start, params_norm.min_N₀, params_norm.max_N₀)
    println("[start] x₀ = $(x₀) (clamped to [$(params_norm.min_N₀), $(params_norm.max_N₀)])")
    f₀ = obj(x₀)
    println("[start] χ²(x₀) = $(f₀)")

    xᵣ, fᵣ = x₀, f₀
    n = 0
    while n < params_norm.max_expansions
        xtry = clamp(xᵣ * params_norm.step_factor, params_norm.min_N₀, params_norm.max_N₀)
        if xtry == xᵣ; println("[up] no further expansion (xtry == xr = $xᵣ)"); break; end
        println("[up] trying x = $(xtry) (factor=$(params_norm.step_factor))")
        ftry = obj(xtry)
        println("[up] χ²($(xtry)) = $(ftry)")
        if ftry < fᵣ
            println("[up] improved: $(fᵣ) → $(ftry); updating xr")
            xᵣ, fᵣ = xtry, ftry
            n += 1
        else
            println("[up] no improvement; stop upward expansion")
            break
        end
    end

    xₗ, fₗ = x₀, f₀
    n = 0
    while n < params_norm.max_expansions
        xtry = clamp(xₗ / params_norm.step_factor, params_norm.min_N₀, params_norm.max_N₀)
        if xtry == xₗ; println("[down] no further expansion (xtry == xl = $xₗ)"); break; end
        println("[down] trying x = $(xtry) (factor=1/$(params_norm.step_factor))")
        ftry = obj(xtry)
        println("[down] χ²($(xtry)) = $(ftry)")
        if ftry < fₗ
            println("[down] improved: $(fₗ) → $(ftry); updating xl")
            xₗ, fₗ = xtry, ftry
            n += 1
        else
            println("[down] no improvement; stop downward expansion")
            break
        end
    end

    a = min(xₗ, xᵣ); b = max(xₗ, xᵣ)
    if a == b
        a = max(params_norm.min_N₀, x₀/params_norm.step_factor)
        b = min(params_norm.max_N₀, x₀*params_norm.step_factor)
        println("[bracket] a==b fallback → a=$(a), b=$(b)")
    end
    println("[bracket] optimizing in interval: a = $(a), b = $(b)")
    if a >= b
        error("x_lower must be less than x_upper: a = $a, b = $b")
    end

    println("[cache] evaluations so far: $(length(cache)) entries")
    for (k,v) in sort(collect(cache); by = x->x[1])
        println("  cache: N₀=$(k) → χ²=$(v)")
    end

    result = optimize(obj, a, b, Brent(); rel_tol=params_norm.brent_reltol)
    best_N₀ = Optim.minimizer(result)
    best_χsq = Optim.minimum(result)

    println("[result] best N₀ = $(best_N₀) with χ² = $(best_χsq)")

     return best_N₀
end

function find_best_N₀_at_Δ₀(params_shape, params_wavefct, params_mc, params_run, params_norm)

    lo = params_norm.min_N₀
    hi = params_norm.max_N₀
    factor = get(params_norm, :step_factor, 1.5)
    start = clamp(get(params_norm, :start, 1.0), lo, hi)
    max_steps = get(params_norm, :max_expansions, 10)
    ngrid = get(params_norm, :ngrid, 21)

    cache = Dict{Float64,Float64}()

    obj(N) = get!(cache, N) do
        println("[eval] χ² at N₀ = $(N)")
        χsq_for_N₀_at_Δ₀(N, params_shape, params_wavefct, params_mc, params_run, params_norm)
    end

    s = clamp(start, lo, hi)
    println("[start] start = $(s) (clamped to [$(lo), $(hi)])")
    f_s = obj(s)
    println("[start] χ²(start) = $(f_s)")

    best_up = s
    f_up = f_s
    x = s
    for i in 1:max_steps
        xtry = clamp(x * factor, lo, hi)
        if xtry == x
            println("[up] reached boundary at x=$(x); stop")
            break
        end
        println("[up] trying x = $(xtry)")
        ftry = obj(xtry)
        println("[up] χ²($(xtry)) = $(ftry)")
        if ftry < f_up
            best_up, f_up = xtry, ftry
            x = xtry
        else
            println("[up] no improvement at x=$(xtry); stop upward expansion")
            break
        end
    end

    best_down = s
    f_down = f_s
    x = s
    for i in 1:max_steps
        xtry = clamp(x / factor, lo, hi)
        if xtry == x
            println("[down] reached boundary at x=$(x); stop")
            break
        end
        println("[down] trying x = $(xtry)")
        ftry = obj(xtry)
        println("[down] χ²($(xtry)) = $(ftry)")
        if ftry < f_down
            best_down, f_down = xtry, ftry
            x = xtry
        else
            println("[down] no improvement at x=$(xtry); stop downward expansion")
            break
        end
    end

    a = max(lo, min(best_down, best_up, s))
    b = min(hi, max(best_down, best_up, s))
    if a == b
        a = max(lo, s / factor)
        b = min(hi, s * factor)
        println("[bracket] a==b fallback → a=$(a), b=$(b)")
    end

    println("[bracket] scanning linear grid between a=$(a) and b=$(b) with ngrid=$(ngrid)")
    N₀_values = collect(range(a, stop=b, length=ngrid))

    if all(abs.(N₀_values .- s) .> eps(Float64))
        push!(N₀_values, s)
        sort!(N₀_values)
        N₀_values = unique(N₀_values)
    end

    N₀_list = Float64[]
    χsq_list = Float64[]
    best_χsq = Inf
    best_N₀ = nothing

    for N₀_val in N₀_values
        f = try
            obj(N₀_val)
        catch e
            @warn "evaluation failed for N₀ = $N₀_val; setting χ² = Inf" exception=(e,)
            Inf
        end

        push!(N₀_list, N₀_val)
        push!(χsq_list, f)

        if isfinite(f) && f < best_χsq
            best_χsq = f
            best_N₀ = N₀_val
        end
    end

    println("[result] best N₀ = $(best_N₀) with χ² = $(best_χsq)")
    return best_N₀, best_χsq, N₀_list, χsq_list
end

function χsq_from_saved(N₀, params_run, params_norm)
    try
        p_run_iter = if params_norm.unique_outdirs
            outdir = joinpath(params_run.outdir, "N0_$(round(N₀; digits=6))")
            merge(params_run, (outdir = outdir,))
        else
            params_run
        end

        compute_dir = joinpath(p_run_iter.savepath, p_run_iter.outdir)
        endswith(compute_dir, "/") && (compute_dir = compute_dir[1:end-1])

        Nsamples = params_norm.nsamples_norm
        Δ0 = params_norm.Δ₀

        # read precomputed cross sections from compute_dir
        t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err =
            compute_cross_sections(compute_dir, [Δ0], Nsamples, p_run_iter, Nsamples)

        tcent_coh_hera, dσdt_coh_hera, Δtot_coh_hera = read_coherent_data(params_norm.coherent_data_path)

        model_val = dσdt_coh[1]
        data_val  = dσdt_coh_hera[1]
        data_err  = Δtot_coh_hera[1]

        χsq = ((model_val - data_val) / data_err)^2
        return isfinite(χsq) ? χsq : Inf
    catch e
        @warn "χsq_from_saved failed" N₀=N₀ exception=e
        return Inf
    end
end

function find_best_N₀_from_saved(N₀_values::AbstractVector{<:Real}, params_run, params_norm)
    best_χsq = Inf
    best_N₀ = nothing
    Nlist = Float64[]
    χlist = Float64[]

    for N0 in N₀_values
        f = χsq_from_saved(N0, params_run, params_norm)
        push!(Nlist, float(N0))
        push!(χlist, f)
        if isfinite(f) && f < best_χsq
            best_χsq = f
            best_N₀ = float(N0)
        end
    end

    return best_N₀, best_χsq, Nlist, χlist
end