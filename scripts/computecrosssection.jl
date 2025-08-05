include("../src/ProtonShape.jl")
using .ProtonShape

using JLD2

function main(params_path::String, nconfigs::Int)
    params_file = joinpath(params_path, "params_used.jld2")
    params = JLD2.load(params_file)
    params_mc = params["p_mc"]
    params_shape = params["p_shape"]
    params_run = params["p_run"]

    Δ_range = range(params_mc.Δmin, stop=params_mc.Δmax, length=params_mc.Δlen)
    Nsamples = params_shape.Nsamples

    t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err = compute_cross_sections(params_path, Δ_range, Nsamples, params_run, nconfigs)

    results_file = joinpath(params_path, "cross_sections.jld2")
    @save results_file t_range dσdt_coh dσdt_coh_err dσdt_incoh dσdt_incoh_err
end

if abspath(PROGRAM_FILE) == @__FILE__
    params_path = length(ARGS) > 0 ? ARGS[1] : error("Usage: julia computecrosssection.jl <path>")
    nconfigs = length(ARGS) > 0 ? parse(Int, ARGS[2]) : 1
    main(params_path, nconfigs)
end