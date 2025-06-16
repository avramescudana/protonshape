using Serialization

include("src/ProtonShape.jl")
using .ProtonShape

t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err = diffractive("coh+incoh", "shapeamp", params_wavefct, params_mc; p_shape=params_shape, run_threads=true)

# Get output file from command line argument, or use default
output_file = length(ARGS) > 0 ? ARGS[1] : "results/test.jls"
save_data = (; t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err,
    params_wavefct, params_gbw, params_cq, params_mc)
serialize(output_file, save_data)
