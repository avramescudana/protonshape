# include("../src/ProtonShape.jl")
# using .ProtonShape

import Pkg

function find_project_root(startdir = @__DIR__)
    dir = normpath(startdir)
    while true
        if isfile(joinpath(dir, "Project.toml")) || isfile(joinpath(dir, "src", "ProtonShape.jl"))
            return dir
        end
        parent = dirname(dir)
        if parent == dir
            error("Could not find project root containing Project.toml or src/ProtonShape.jl")
        end
        dir = parent
    end
end

proj_root = find_project_root()
Pkg.activate(proj_root)

# include(joinpath(proj_root, "src", "ProtonShape.jl"))
# using .ProtonShape
using ProtonShape

using Serialization

using JLD2

function main(params_path::String, nconfigs::Int, multiple_configs::Int)
    if multiple_configs == 1
        subdirs = filter(f -> isdir(joinpath(params_path, f)), readdir(params_path))
        isempty(subdirs) && error("No subdirectories found in $(params_path) with multiple_configs == 1")

        println("Subdirectories: $subdirs")

        for subdir in sort(subdirs)
            println("Processing subdirectory: $subdir")
            subdir_path = joinpath(params_path, subdir)
            params_file = joinpath(subdir_path, "params.jld2")

            params = JLD2.load(params_file)
            params_mc = params["params_mc"]
            params_shape = params["params_shape_eff_best_N₀"]
            params_run = params["params_run_sigma"]
            # loaded = open(params_file, "r") do io
            #     deserialize(io)
            # end

            # (diff, dip, params_wavefc2, params_mc, params_gbw, params_cq, params_shape, params_run) = loaded

            # Δ_range = range(params_mc.Δmin, stop=params_mc.Δmax, length=params_mc.Δlen)
            Δ_range = collect(range(params_mc.Δmin, stop=params_mc.Δmax, length=params_mc.Δlen))
            Nsamples = params_shape.Nsamples

            # t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err = compute_cross_sections(
            #     subdir_path, Δ_range, Nsamples, params_run, nconfigs)
            t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err =
                compute_cross_sections(subdir_path, Δ_range, Nsamples, params_run, nconfigs)

            results_file = joinpath(subdir_path, "cross_sections.jld2")
            @save results_file t_range dσdt_coh dσdt_coh_err dσdt_incoh dσdt_incoh_err
        end
    else
        params_file = joinpath(params_path, "params.jld2")
        params = JLD2.load(params_file)
        # params_mc = params["p_mc"]
        # params_shape = params["p_shape"]
        # params_run = params["p_run"]
        params_mc = params["params_mc"]
        params_shape = params["params_shape_eff_best_N₀"]
        params_run = params["params_run_sigma"]
        # loaded = open(params_file, "r") do io
        #     deserialize(io)
        # end

        # (diff, dip, params_wavefc2, params_mc, params_gbw, params_cq, params_shape, params_run) = loaded

        # Δ_range = range(params_mc.Δmin, stop=params_mc.Δmax, length=params_mc.Δlen)
        Δ_range = collect(range(params_mc.Δmin, stop=params_mc.Δmax, length=params_mc.Δlen))
        Nsamples = params_shape.Nsamples

        # t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err = compute_cross_sections(
        #     params_path, Δ_range, Nsamples, params_run, nconfigs, multiple_configs)
        t_range, dσdt_coh, dσdt_coh_err, dσdt_incoh, dσdt_incoh_err =
            compute_cross_sections(params_path, Δ_range, Nsamples, params_run, nconfigs)

        results_file = joinpath(params_path, "cross_sections.jld2")
        @save results_file t_range dσdt_coh dσdt_coh_err dσdt_incoh dσdt_incoh_err
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    params_path = length(ARGS) > 0 ? ARGS[1] : error("Usage: julia computecrosssection.jl <path> <nconfigs> <multiple_configs>")
    nconfigs = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 1
    multiple_configs = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 1
    println("Running main with: $params_path, $nconfigs configs, multiple_configs=$multiple_configs")
    main(params_path, nconfigs, multiple_configs)
end
