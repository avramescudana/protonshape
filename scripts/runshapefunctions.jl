using Serialization
using Random, UUIDs

include("../src/ProtonShape.jl")
using .ProtonShape

arrayindex = parse(Int, ARGS[1])
nconfigs   = parse(Int, ARGS[2])
m          = parse(Int, ARGS[3])
nmax       = parse(Int, ARGS[4])
savepath   = ARGS[5]
sigma_list = [parse(Float64, s) for s in ARGS[6:end]]

println("arrayindex = ", arrayindex)
println("nconfigs   = ", nconfigs)
println("m          = ", m)
println("nmax       = ", nmax)
println("savepath   = ", savepath)
println("sigma_list = ", sigma_list)

Random.seed!(arrayindex)

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
    diffractive("coh+incoh", "shapeamp", params_wavefct, params_mc; 
                p_shape=params_shape_sigma, p_run=params_run_sigma)
end
