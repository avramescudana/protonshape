using Serialization
using Random, UUIDs

include("../src/ProtonShape.jl")
using .ProtonShape

arrayindex = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 1
nconfigs = length(ARGS) > 0 ? parse(Int, ARGS[2]) : 1
params_run_array = merge(params_run, (arrayindex = arrayindex,))

params_shape_eff = merge(params_shape, (Nsamples = nconfigs,))

if params_shape.type=="samemn"
    coeff_dicts = sample_amp_dict_same_mn(params_shape_eff)
elseif params_shape.type=="samem_multin"
    coeff_dicts = sample_amp_dict_samem_multin(params_shape_eff)
else
    error("Unknown sampling type: $(p_shape.type)")
end

params_run_eff = merge(params_run_array, (amp_dict = coeff_dicts,))

diffractive("coh+incoh", "shapeamp", params_wavefct, params_mc; p_shape=params_shape_eff, p_run=params_run_eff)
