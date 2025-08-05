using Serialization
using Random, UUIDs

include("../src/ProtonShape.jl")
using .ProtonShape

arrayindex = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 1
nconfigs = length(ARGS) > 0 ? parse(Int, ARGS[2]) : 1
sigma_idx  = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 1

σ_values = range(1.0, stop=20.0, length=10)
σ_val = σ_values[sigma_idx]

params_shape_eff = merge(params_shape, (Nsamples = nconfigs, σ = σ_val,))
outdir_name = "sigma_$(σ_val)"
params_run_array = merge(params_run, (arrayindex = arrayindex, outdir = outdir_name))

if params_shape.type=="samemn"
    coeff_dicts = sample_amp_dict_same_mn(params_shape_eff)
elseif params_shape.type=="samem_multin"
    coeff_dicts = sample_amp_dict_samem_multin(params_shape_eff)
else
    error("Unknown sampling type: $(p_shape.type)")
end

params_run_eff = merge(params_run_array, (amp_dict = coeff_dicts,))

diffractive("coh+incoh", "shapeamp", params_wavefct, params_mc; p_shape=params_shape_eff, p_run=params_run_eff)
