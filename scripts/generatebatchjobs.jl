using Printf
using Serialization
using Random

if length(ARGS) < 1
    error("Usage: julia generatebatchjobs.jl path/to/paramsfile.jl [mode]\n  mode: csc (default) | local")
end

params_file = ARGS[1]
mode = length(ARGS) >= 2 ? ARGS[2] : "csc"
mode ∉ ("csc", "local") && error("Unknown mode: $mode — must be 'csc' or 'local'")

params_basename = split(basename(params_file), ".")[1]
output_dir = joinpath("scripts", "sbatchjobs")
mkpath(output_dir)

output_file =
    mode == "csc" ? joinpath(output_dir, "submitalljobs_$(params_basename).sh") :
                    joinpath(output_dir, "runalljobs_$(params_basename)_local.sh")

include(params_file)
using .SimulationParams

open(output_file, "w") do io
    write(io, "#!/bin/bash\n\n")
    write(io, "set -euo pipefail\n\n")

    for (set_index, params) in enumerate(SimulationParams.param_sets)
        m          = params.m
        nmax       = params.nmax
        sigma_list = params.sigma_list
        savepath   = params.savepath
        nconfigs   = params.nconfigs
        N₀_list    = params.N₀_list
        paramset   = params.paramset
        find_norm  = params.find_norm
        nsamples_norm = get(params, :nsamples_norm, 20)

        overrides = get(params, :overrides, nothing)
        if overrides !== nothing
            mkpath(savepath) 
            override_file = joinpath(savepath, "params_override.jl2")
            open(override_file, "w") do f
                write(f, "const params_override = $(repr(overrides))\n")
            end
        end

        for sigma in sigma_list
            if mode == "csc"
                if find_norm
                    norm_randomseed = rand(Int64)
                    sbatch_norm_command =
"""sbatch <<'EOF'
#!/bin/bash
#SBATCH --job-name=runshapefunctions_norm_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=$(savepath)/slurm_out/runshapefunctions_norm_$(set_index)_$(paramset)_$(sigma).out
#SBATCH --error=$(savepath)/slurm_out/runshapefunctions_norm_$(set_index)_$(paramset)_$(sigma).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4000

module load julia

julia --project=. scripts/runshapefunctions.jl \\
    $(1) \\
    $(nsamples_norm) \\
    $(norm_randomseed) \\
    $(m) \\
    $(nmax) \\
    $(savepath) \\
    $(sigma) \\
    true \\
    -1.0 \\
    $(paramset) \\
    norm_only
EOF

"""
                    write(io, sbatch_norm_command * "\n")
                end

                for N₀ in N₀_list
                    for config_index in 1:nconfigs
                        randomseed = rand(Int64)
                        placeholder_N₀ = find_norm ? -1.0 : N₀

                        sbatch_command =
"""sbatch <<'EOF'
#!/bin/bash
#SBATCH --job-name=runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index)
#SBATCH --output=$(savepath)/slurm_out/runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index).out
#SBATCH --error=$(savepath)/slurm_out/runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000

module load julia

julia --project=. scripts/runshapefunctions.jl \\
    $(config_index) \\
    $(nconfigs) \\
    $(randomseed) \\
    $(m) \\
    $(nmax) \\
    $(savepath) \\
    $(sigma) \\
    false \\
    $(placeholder_N₀) \\
    $(paramset) \\
    run
EOF

"""
                        write(io, sbatch_command)
                    end
                end
            else
                if find_norm
                    norm_randomseed = rand(Int64)
                    cmd = @sprintf("echo 'Running norm_only: set=%d sigma=%s' \nmkdir -p %s\njulia --project=. scripts/runshapefunctions.jl %d %d %d %d %d %s %g true %g %s norm_only\n\n",
                        set_index, string(sigma), savepath, 1, nsamples_norm, norm_randomseed, m, nmax, savepath, sigma, -1.0, paramset)
                    write(io, cmd)
                end

                for N₀ in N₀_list
                    for config_index in 1:nconfigs
                        randomseed = rand(Int64)
                        placeholder_N₀ = find_norm ? -1.0 : N₀
                        cmd = @sprintf("echo 'Running job: set=%d paramset=%s sigma=%g config=%d' \nmkdir -p %s\njulia --project=. scripts/runshapefunctions.jl %d %d %d %d %d %s %g false %g %s run\n\n",
                                set_index, paramset, sigma, config_index, savepath, config_index, nconfigs, randomseed, m, nmax, savepath, sigma, placeholder_N₀, paramset)
                        write(io, cmd)
                    end
                end
            end
        end
    end
end

println("Generated $mode jobs script: $output_file")