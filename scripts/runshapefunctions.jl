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

println("arrayindex = ", arrayindex)
println("nconfigs   = ", nconfigs)
println("randomseed = ", randomseed)
println("m          = ", m)
println("nmax       = ", nmax)
println("savepath   = ", savepath)
println("sigma      = ", sigma)
println("find_norm  = ", find_norm)
println("N₀         = ", N₀)

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
    params_shape_eff_best_N₀ = merge(params_shape_eff, (N₀ = N₀,))
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