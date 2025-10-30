using Serialization
using Random, UUIDs
using JLD2

import Pkg
Pkg.activate("/users/davrames/protonshape")

include("/users/davrames/protonshape/src/ProtonShape.jl")
using .ProtonShape

arrayindex = parse(Int, ARGS[1])
nconfigs   = parse(Int, ARGS[2])
randomseed = parse(Int64, ARGS[3])  
m          = parse(Int, ARGS[4])
nmax       = parse(Int, ARGS[5])
savepath   = ARGS[6]
sigma      = parse(Float64, ARGS[7])
find_norm  = parse(Bool, ARGS[8])
N₀         = parse(Float64, ARGS[9])
paramset   = length(ARGS) >= 10 ? ARGS[10] : "default_paramset"
mode       = length(ARGS) >= 11 ? ARGS[11] : "run"

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

# Deterministic location for norm file
norm_dir = joinpath(savepath, "norm")
mkpath(norm_dir)
norm_file = joinpath(norm_dir, string(paramset, "_sigma_", sigma, "_bestN0.jld2"))

if mode == "norm_only"
    # Run the adaptive search for best N₀ using params_shape_eff (Nsamples = nsamples_norm passed in as nconfigs)
    best_N₀ = find_best_N₀_at_Δ₀_adaptive(
        params_shape_eff,
        params_wavefct,
        params_mc,
        params_run_eff_base,   # use base run params while searching
        params_norm
    )
    println("find_norm -> best_N₀ = ", best_N₀)

    # Save best N₀
    @save norm_file best_N₀ params_shape_eff
    println("Saved best N₀ to $norm_file")
    # Exit without doing the full diffractive run
    exit()
end

# mode == "run"
# If N₀ is a placeholder (<= 0), wait for norm_file and load best_N₀
if N₀ <= 0
    println("No N₀ provided; attempting to load best_N₀ from norm file: $norm_file")
    max_tries = 3*360    # 360 * 10s = 3600s = 1 hour
    found = false
    for i in 1:max_tries
        if isfile(norm_file)
            println("Found norm file: $norm_file")
            data = JLD2.load(norm_file)
            if haskey(data, "best_N₀")
                best_N₀ = data["best_N₀"]
                println("Loaded best_N₀ = ", best_N₀)
                params_shape_eff = merge(params_shape_eff, (N₀ = best_N₀,))
                found = true
                break
            else
                println("norm file exists but does not contain best_N₀, retrying...")
            end
        end
        sleep(10)  # wait 10s before retry
    end
    if !found
        error("Timed out waiting for normalization file: $norm_file")
    end
else
    # Use provided N₀
    params_shape_eff = merge(params_shape_eff, (N₀ = N₀,))
end

# Prepare coeff_dicts as before
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
    # This branch remains for backward compatibility: if find_norm is true in a run job, compute & use best N₀ locally.
    best_N₀ = find_best_N₀_at_Δ₀_adaptive(
        params_shape_eff,
        params_wavefct,
        params_mc,
        params_run_eff_base,   # use base run params while searching
        params_norm
    )
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