using Serialization
using Random, UUIDs
using JLD2

include("../src/ProtonShape.jl")
using .ProtonShape

arrayindex = parse(Int, ARGS[1])
nconfigs   = parse(Int, ARGS[2])
randomseed = parse(Int64, ARGS[3])  
m          = parse(Int, ARGS[4])
nmax       = parse(Int, ARGS[5])
savepath   = ARGS[6]
sigma      = parse(Float64, ARGS[7])
N₀         = parse(Float64, ARGS[8])

println("arrayindex = ", arrayindex)
println("nconfigs   = ", nconfigs)
println("randomseed = ", randomseed)
println("m          = ", m)
println("nmax       = ", nmax)
println("savepath   = ", savepath)
println("sigma      = ", sigma)
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

# The rest of your logic, now for a single sigma and N₀:
coeff_dicts =
    if params_shape_eff.type == "samemn"
        sample_amp_dict_same_mn(params_shape_eff)
    elseif params_shape_eff.type == "samem_multin"
        sample_amp_dict_samem_multin(params_shape_eff)
    else
        error("Unknown sampling type: $(params_shape_eff.type)")
    end

outdir_name = "sigma_$(sigma)_N0_$(N₀)"
params_run_sigma = merge(params_run_eff_base, (
    outdir = outdir_name,
    amp_dict = coeff_dicts,
))

diff = "coh+incoh"
dip = "shapeamp"
diffractive(diff, dip, params_wavefct, params_mc; p_shape=params_shape_eff, p_run=params_run_sigma)

if arrayindex == 1
    outdir_path = params_run_sigma.savepath * "/" * params_run_sigma.outdir
    params_file = joinpath(outdir_path, "params.jld2")
    @save params_file diff dip params_wavefct params_mc params_shape_eff params_run_sigma
end