using Serialization
using Random

include("../src/ProtonShape.jl")
using .ProtonShape

# t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err = diffractive("coh+incoh", "shapeamp", params_wavefct, params_mc; p_shape=params_shape, run_threads=false)

arrayindex = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 1
params_run_updated = merge(params_run, (arrayindex = arrayindex,))
Random.seed!(12345 + arrayindex)

diffractive("coh+incoh", "shapeamp", params_wavefct, params_mc; p_shape=params_shape, p_run=params_run_updated)

# Get output file from command line argument, or use default
# output_file = length(ARGS) > 0 ? ARGS[1] : "results/test.jls"
# save_data = (; t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err,
    # params_wavefct, params_gbw, params_cq, params_mc)
# serialize(output_file, save_data)
