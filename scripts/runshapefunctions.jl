using Serialization
using Random, UUIDs
using JLD2

import Pkg

function find_project_root(startdir = @__DIR__)
    dir = normpath(startdir)
    while true
        if isfile(joinpath(dir, "Project.toml")) || isfile(joinpath(dir, "src", "ProtonShape.jl"))
            return dir
        end
        parent = dirname(dir)
        if parent == dir
            error("Could not find project root containing Project.toml or src/ProtonShape.jl")
        end
        dir = parent
    end
end

proj_root = find_project_root()
Pkg.activate(proj_root)

include(joinpath(proj_root, "src", "ProtonShape.jl"))
using .ProtonShape

# ----------------- parse arguments -----------------

arrayindex = parse(Int, ARGS[1])
arg2       = parse(Int, ARGS[2]) # nconfigs OR nsamples_norm depending on mode
randomseed = parse(Int64, ARGS[3])
m          = parse(Int, ARGS[4])
nmax       = parse(Int, ARGS[5])
savepath   = ARGS[6]
sigma      = parse(Float64, ARGS[7])
find_norm  = parse(Bool, ARGS[8])
N₀         = parse(Float64, ARGS[9])
paramset   = length(ARGS) >= 10 ? ARGS[10] : "default_paramset"
mode       = length(ARGS) >= 11 ? ARGS[11] : "run"

# Determine nconfigs vs nsamples_norm consistently
if mode == "norm_only" || mode == "norm_point"
    nsamples_norm = arg2
    nconfigs      = nsamples_norm
else
    nconfigs      = arg2
    nsamples_norm = get(params_norm, :nsamples_norm, 1)
end

# ----------------- apply overrides if present -----------------

override_file = joinpath(savepath, "params_override.jl2")

if isfile(override_file)
    include(override_file)
    if isdefined(Main, :params_override)
        for (k, v) in pairs(params_override)
            if isdefined(ProtonShape, k)
                base   = getfield(ProtonShape, k)
                merged = merge(base, v)
                @info "Applying override for $(k)" v
                @eval Main $(k) = $merged
            else
                @warn "Unknown override key: $(k) (skipping)"
            end
        end
    end
end

println("arrayindex = ", arrayindex)
println("nconfigs   = ", nconfigs)
println("randomseed = ", randomseed)
println("m          = ", m)
println("nmax       = ", nmax)
println("savepath   = ", savepath)
println("sigma      = ", sigma)
println("find_norm  = ", find_norm)
println("N₀         = ", N₀)
println("paramset   = ", paramset)
println("mode       = ", mode)

Random.seed!(randomseed)

params_shape_eff = merge(params_shape, (
    Nsamples = nconfigs,
    mn       = (m, 1),
    nvals    = nmax,
    N₀       = N₀,
    σ        = sigma,
))

params_run_eff_base = merge(params_run, (
    savepath   = savepath,
    arrayindex = arrayindex,
))

# For normalization runs, base outdir = "norm"
params_run_eff_base_norm = merge(params_run_eff_base, (
    outdir = "norm",
))

norm_dir  = joinpath(savepath, "norm")
mkpath(norm_dir)
norm_file = joinpath(norm_dir, string(paramset, "_sigma_", sigma, "_bestN0.jld2"))

# choose params_mc for this run
params_mc_eff =
    if mode == "norm_only" || mode == "norm_point"
        Δ₀ = params_norm.Δ₀
        merge(params_mc, (Δmin = Δ₀, Δmax = Δ₀, Δlen = 1))
    else
        params_mc
    end

# ----------------- norm_only -----------------

if mode == "norm_only" && isfile(norm_file)
    data = JLD2.load(norm_file)
    if haskey(data, "best_N₀")
        println("Norm file already exists; using existing best_N₀ = ", data["best_N₀"], " — skipping normalization")
        exit()
    else
        println("Norm file exists but missing best_N₀ — will recompute")
    end
end

if mode == "norm_only"
    params_norm_eff = merge(params_norm, (nsamples_norm = nsamples_norm,))
    best_N₀, best_χsq, _, _ =
        find_best_N₀_at_Δ₀(params_shape_eff, params_wavefct, params_mc_eff, params_run_eff_base_norm, params_norm_eff)
    println("find_norm (norm_only) -> best_N₀ = ", best_N₀)

    @save norm_file best_N₀ params_shape_eff
    println("Saved best N₀ to $norm_file")
    exit()
end

# new: write a deterministic N0 grid (no MC) and save to plain text for bash
if mode == "norm_grid"
    println("Generating N₀ grid (norm_grid mode)")

    lo = params_norm.min_N₀
    hi = params_norm.max_N₀
    start = params_norm.start
    ngrid = get(params_norm, :ngrid, 21)
    step_factor = get(params_norm, :step_factor, 1.5)

    # simple bracket around start (deterministic, no χ² evals)
    a = max(lo, start / step_factor)
    b = min(hi, start * step_factor)
    if a >= b
        # fallback to full linear range
        a, b = lo, hi
    end

    N0_values = collect(range(a, stop=b, length=ngrid))
    if all(abs.(N0_values .- start) .> eps(Float64))
        push!(N0_values, start)
        sort!(N0_values)
        N0_values = unique(N0_values)
    end

    mkpath(norm_dir)
    grid_jld = joinpath(norm_dir, string(paramset, "_sigma_", sigma, "_N0grid.jld2"))
    grid_txt = joinpath(norm_dir, string(paramset, "_sigma_", sigma, "_N0grid.txt"))

    @save grid_jld N0_values
    open(grid_txt, "w") do f
        write(f, join(string.(N0_values), " "))
    end

    println("Wrote N₀ grid files: ", grid_jld, "  ", grid_txt)
    exit()
end

# modify norm_collect to accept a grid txt path (ARGS[12]) or default to naming convention
if mode == "norm_collect"
    # ARGS[12] (optional) = path to N0 grid txt file
    grid_txt = length(ARGS) >= 12 ? ARGS[12] : joinpath(norm_dir, string(paramset, "_sigma_", sigma, "_N0grid.txt"))
    if !isfile(grid_txt)
        error("norm_collect: grid file not found: $grid_txt")
    end
    N0_values = parse.(Float64, split(chomp(read(grid_txt, String))))
    println("norm_collect: reading grid from $grid_txt -> $(length(N0_values)) entries")

    # use existing ProtonShape helpers (we added find_best_N₀_from_saved earlier)
    best_N₀, best_χsq, _, _ = ProtonShape.find_best_N₀_from_saved(N0_values, params_run_eff_base_norm, params_norm)

    if best_N₀ === nothing
        error("norm_collect: failed to determine best_N₀")
    end

    println("norm_collect -> best_N₀ = ", best_N₀, " χ² = ", best_χsq)
    @save norm_file best_N₀ best_χsq
    println("Saved best N₀ to $norm_file")
    exit()
end

# ----------------- norm_grid -----------------

if mode == "norm_grid"
    println("Generating N₀ grid (norm_grid mode)")

    lo          = params_norm.min_N₀
    hi          = params_norm.max_N₀
    start       = params_norm.start
    ngrid       = get(params_norm, :ngrid, 21)
    step_factor = get(params_norm, :step_factor, 1.5)

    a = max(lo, start / step_factor)
    b = min(hi, start * step_factor)
    if a >= b
        a, b = lo, hi
    end

    N0_values = collect(range(a, stop=b, length=ngrid))
    if all(abs.(N0_values .- start) .> eps(Float64))
        push!(N0_values, start)
        sort!(N0_values)
        N0_values = unique(N0_values)
    end

    mkpath(norm_dir)
    grid_jld = joinpath(norm_dir, string(paramset, "_sigma_", sigma, "_N0grid.jld2"))
    grid_txt = joinpath(norm_dir, string(paramset, "_sigma_", sigma, "_N0grid.txt"))

    @save grid_jld N0_values
    open(grid_txt, "w") do f
        write(f, join(string.(N0_values), " "))
    end

    println("Wrote N₀ grid files: ", grid_jld, "  ", grid_txt)
    exit()
end

# ----------------- norm_collect -----------------

if mode == "norm_collect"
    grid_txt = length(ARGS) >= 12 ?
        ARGS[12] :
        joinpath(norm_dir, string(paramset, "_sigma_", sigma, "_N0grid.txt"))

    if !isfile(grid_txt)
        error("norm_collect: grid file not found: $grid_txt")
    end
    N0_values = parse.(Float64, split(chomp(read(grid_txt, String))))
    println("norm_collect: reading grid from $grid_txt -> $(length(N0_values)) entries")

    best_N₀, best_χsq, _, _ = ProtonShape.find_best_N₀_from_saved(N0_values, params_run_eff_base_norm, params_norm)

    if best_N₀ === nothing
        error("norm_collect: failed to determine best_N₀")
    end

    println("norm_collect -> best_N₀ = ", best_N₀, " χ² = ", best_χsq)
    @save norm_file best_N₀ best_χsq
    println("Saved best N₀ to $norm_file")
    exit()
end

# ----------------- core pipeline (norm_point and run) -----------------

# For run: if N₀ <= 0 we wait for norm_file and load best_N₀
if mode == "run" && N₀ <= 0
    println("No N₀ provided; attempting to load best_N₀ from norm file: $norm_file")
    max_tries = 3 * 360
    global found = false
    for i in 1:max_tries
        if isfile(norm_file)
            println("Found norm file: $norm_file")
            data = JLD2.load(norm_file)
            if haskey(data, "best_N₀")
                best_N₀ = data["best_N₀"]
                println("Loaded best_N₀ = ", best_N₀)
                global params_shape_eff = merge(params_shape_eff, (N₀ = best_N₀,))
                global find_norm = false
                global params_shape_eff_best_N₀ = params_shape_eff
                global found = true
                break
            else
                println("norm file exists but does not contain best_N₀, retrying...")
            end
        end
        sleep(10)
    end
    if !found
        error("Timed out waiting for normalization file: $norm_file")
    end
elseif mode == "run"
    global params_shape_eff = merge(params_shape_eff, (N₀ = N₀,))
end

# For norm_point, N₀ is passed directly and we don't touch it here.

# Sample amplitudes
coeff_dicts =
    if params_shape_eff.type == "samemn"
        sample_amp_dict_same_mn(params_shape_eff)
    elseif params_shape_eff.type == "samem_multin"
        sample_amp_dict_samem_multin(params_shape_eff)
    else
        error("Unknown sampling type: $(params_shape_eff.type)")
    end

params_shape_eff_best_N₀ = params_shape_eff

diff = "coh+incoh"
dip  = "shapeamp"

# ---------- norm_point: multiple configs per N₀ ----------

if mode == "norm_point"
    base_outdir = params_run_eff_base_norm.outdir   # "norm"
    N0_tag      = "N0_$(round(N₀; digits=6))"
    full_outdir = joinpath(base_outdir, N0_tag)

    params_run_sigma_base = merge(params_run_eff_base_norm, (
        outdir   = full_outdir,
        amp_dict = coeff_dicts,
    ))

    println("norm_point: N₀ = $N₀, nsamples_norm = $nsamples_norm → generating configs 1:$nsamples_norm")

    for cfg in 1:nsamples_norm
        p_run_cfg = merge(params_run_sigma_base, (arrayindex = cfg,))
        println("norm_point: running cfg = $cfg in outdir = $(p_run_cfg.outdir)")
        diffractive(diff, dip, params_wavefct, params_mc_eff; p_shape=params_shape_eff_best_N₀, p_run=p_run_cfg)

        if cfg == 1
            outdir_path = p_run_cfg.savepath * "/" * p_run_cfg.outdir
            mkpath(outdir_path)
            params_file = joinpath(outdir_path, "params.jld2")
            @save params_file diff dip params_wavefct params_mc params_shape_eff_best_N₀ p_run_cfg
        end
    end

    exit()
end

# ---------- Normal physics run (mode == "run") ----------

params_run_sigma = begin
    outdir_name = "sigma_$(sigma)_N0_$(params_shape_eff_best_N₀.N₀)"
    merge(params_run_eff_base, (
        outdir   = outdir_name,
        amp_dict = coeff_dicts,
    ))
end

diffractive(diff, dip, params_wavefct, params_mc_eff; p_shape=params_shape_eff_best_N₀, p_run=params_run_sigma)

if arrayindex == 1
    outdir_path = params_run_sigma.savepath * "/" * params_run_sigma.outdir
    mkpath(outdir_path)
    params_file = joinpath(outdir_path, "params.jld2")
    @save params_file diff dip params_wavefct params_mc params_shape_eff_best_N₀ params_run_sigma
end
