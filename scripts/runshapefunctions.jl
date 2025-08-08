using Serialization
using Random, UUIDs

include("../src/ProtonShape.jl")
using .ProtonShape

arrayindex = parse(Int, ARGS[1])
nconfigs   = parse(Int, ARGS[2])
sigma_list = [parse(Float64, s) for s in ARGS[3:end]]

Random.seed!(arrayindex)

params_shape_eff = merge(params_shape, (Nsamples = nconfigs,))
if params_shape.type == "samemn"
    coeff_dicts = sample_amp_dict_same_mn(params_shape_eff)
elseif params_shape.type == "samem_multin"
    coeff_dicts = sample_amp_dict_samem_multin(params_shape_eff)
else
    error("Unknown sampling type: $(p_shape.type)")
end

for σ_val in sigma_list
    params_shape_sigma = merge(params_shape_eff, (σ = σ_val,))
    outdir_name = "sigma_$(σ_val)"
    params_run_array = merge(params_run, (arrayindex = arrayindex, outdir = outdir_name))
    params_run_eff = merge(params_run_array, (amp_dict = coeff_dicts,))
    diffractive("coh+incoh", "shapeamp", params_wavefct, params_mc; 
                p_shape=params_shape_sigma, p_run=params_run_eff)
end
