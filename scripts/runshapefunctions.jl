using Serialization
using Random, UUIDs
using JLD2

import Pkg
# Pkg.activate("/users/davrames/protonshape")
# include("/users/davrames/protonshape/src/ProtonShape.jl")

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

arrayindex = parse(Int, ARGS[1])
arg2 = parse(Int, ARGS[2]) # nconfigs OR nsamples_norm depending on mode
randomseed = parse(Int64, ARGS[3])
m          = parse(Int, ARGS[4])
nmax       = parse(Int, ARGS[5])
savepath   = ARGS[6]
sigma      = parse(Float64, ARGS[7])
find_norm  = parse(Bool, ARGS[8])
N₀         = parse(Float64, ARGS[9])
paramset   = length(ARGS) >= 10 ? ARGS[10] : "default_paramset"
mode       = length(ARGS) >= 11 ? ARGS[11] : "run"

# Determine nconfigs vs nsamples_norm consistently:
if mode == "norm_only"
    nsamples_norm = arg2
    nconfigs = nsamples_norm      # use same value when building params_shape_eff for normalization runs
else
    nconfigs = arg2
    nsamples_norm = get(params_norm, :nsamples_norm, 1)  # fallback to default from ProtonShape.params_norm
end

override_file = joinpath(savepath, "params_override.jl2")

if isfile(override_file)
    include(override_file)  
    if isdefined(Main, :params_override)
        for (k, v) in pairs(params_override)   
            if isdefined(ProtonShape, k)
                base = getfield(ProtonShape, k)
                merged = merge(base, v)         
                @info "Applying override for $(k)"
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
    mn = (m, 1),        # initial n will be overwritten below
    nvals = nmax,
    N₀ = N₀,
    σ = sigma,
))

params_run_eff_base = merge(params_run, (
    savepath = savepath,
    arrayindex = arrayindex,
))

params_run_eff_base_norm = merge(params_run_eff_base, (outdir = "sigma_$(sigma)_norm/",))

norm_dir = joinpath(savepath, "norm")
mkpath(norm_dir)
norm_file = joinpath(norm_dir, string(paramset, "_sigma_", sigma, "_bestN0.jld2"))


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
    best_N₀, best_χsq, _, _ = find_best_N₀_at_Δ₀(params_shape_eff, params_wavefct, params_mc, params_run_eff_base_norm, params_norm_eff)
    println("find_norm -> best_N₀ = ", best_N₀)

    @save norm_file best_N₀ params_shape_eff
    println("Saved best N₀ to $norm_file")
    exit()
end

# mode == "run"
# If N₀ is a placeholder (<= 0), wait for norm_file and load best_N₀
if N₀ <= 0
    println("No N₀ provided; attempting to load best_N₀ from norm file: $norm_file")
    max_tries = 3*360    # 360 * 10s = 3600s = 1 hour
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
else
    global params_shape_eff = merge(params_shape_eff, (N₀ = N₀,))
end

coeff_dicts =
    if params_shape_eff.type == "samemn"
        sample_amp_dict_same_mn(params_shape_eff)
    elseif params_shape_eff.type == "samem_multin"
        sample_amp_dict_samem_multin(params_shape_eff)
    else
        error("Unknown sampling type: $(params_shape_eff.type)")
    end

params_shape_eff_best_N₀ = params_shape_eff

if find_norm
    params_norm_eff = merge(params_norm, (nsamples_norm = nsamples_norm,))
    best_N₀ = find_best_N₀_at_Δ₀(params_shape_eff, params_wavefct, params_mc, params_run_eff_base_norm, params_norm_eff)
    println("find_norm -> best_N₀ = ", best_N₀)
    params_shape_eff_best_N₀ = merge(params_shape_eff, (N₀ = best_N₀,))
else
    params_shape_eff_best_N₀ = params_shape_eff
end

outdir_name = "sigma_$(sigma)_N0_$(params_shape_eff_best_N₀.N₀)"
params_run_sigma = merge(params_run_eff_base, (
    outdir = outdir_name,
    amp_dict = coeff_dicts,
))

diff = "coh+incoh"
dip = "shapeamp"

diffractive(diff, dip, params_wavefct, params_mc; p_shape=params_shape_eff_best_N₀, p_run=params_run_sigma)

if arrayindex == 1
    outdir_path = params_run_sigma.savepath * "/" * params_run_sigma.outdir
    params_file = joinpath(outdir_path, "params.jld2")
    @save params_file diff dip params_wavefct params_mc params_shape_eff_best_N₀ params_run_sigma
end