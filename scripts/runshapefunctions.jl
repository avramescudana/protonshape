using Serialization
using Random, UUIDs
using JLD2

include("../src/ProtonShape.jl")
using .ProtonShape

arrayindex = parse(Int, ARGS[1])
nconfigs   = parse(Int, ARGS[2])
randomseed = parse(Int, ARGS[3])  
m          = parse(Int, ARGS[4])
nmax       = parse(Int, ARGS[5])
savepath   = ARGS[6]
sigma_list = [parse(Float64, s) for s in ARGS[7:end]]

println("arrayindex = ", arrayindex)
println("nconfigs   = ", nconfigs)
println("randomseed = ", randomseed)
println("m          = ", m)
println("nmax       = ", nmax)
println("savepath   = ", savepath)
println("sigma_list = ", sigma_list)

Random.seed!(randomseed)

params_shape_eff = merge(params_shape, (
    Nsamples = nconfigs,
    mn = (m, 1),        # initial n will be overwritten below
    nvals = nmax,
))

params_run_eff_base = merge(params_run, (
    savepath = savepath,
    arrayindex = arrayindex,
))

for σ_val in sigma_list
    params_shape_sigma = merge(params_shape_eff, (σ = σ_val,))

    coeff_dicts =
        if params_shape_sigma.type == "samemn"
            sample_amp_dict_same_mn(params_shape_sigma)
        elseif params_shape_sigma.type == "samem_multin"
            sample_amp_dict_samem_multin(params_shape_sigma)
        else
            error("Unknown sampling type: $(params_shape_sigma.type)")
        end

    outdir_name = "sigma_$(σ_val)"
    params_run_sigma = merge(params_run_eff_base, (
        outdir = outdir_name,
        amp_dict = coeff_dicts,
    ))

    diff = "coh+incoh"
    dip = "shapeamp"
    diffractive(diff, dip, params_wavefct, params_mc; p_shape=params_shape_sigma, p_run=params_run_sigma)

    if arrayindex == 1
        outdir_path = params_run_sigma.savepath * "/" * params_run_sigma.outdir
        params_file = joinpath(outdir_path, "params.jld2")
        @save params_file diff dip params_wavefct params_mc params_shape_eff params_run_sigma
    end
end
