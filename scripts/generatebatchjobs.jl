using Printf
using Serialization

if length(ARGS) < 1
    error("Usage: julia generatebatchjobs.jl path/to/paramsfile.jl")
end

params_file = ARGS[1]

params_basename = split(basename(params_file), ".")[1] 
output_file = joinpath("scripts/sbatchjobs", "submitalljobs$(params_basename).sh")

include(params_file)

using .SimulationParams

open(output_file, "w") do io
    write(io, "#!/bin/bash\n\n")

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

        sigma_values_str = join(sigma_list, " ")

        for sigma in sigma_list

            # If this param set requests finding normalization, generate one normalization job per sigma.
            if find_norm
                norm_randomseed = rand(Int64)

                sbatch_norm_command =
                    "sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=runshapefunctions_norm_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=/scratch/lappi/dana/slurm_out/runshapefunctions_norm_$(set_index)_$(paramset)_$(sigma).out
#SBATCH --error=/scratch/lappi/dana/slurm_out/runshapefunctions_norm_$(set_index)_$(paramset)_$(sigma).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4000

module load julia

# mode \"norm_only\" will compute best_N₀ using nsamples_norm and write norm file, then exit.
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
"
                write(io, sbatch_norm_command * "\n")
            end

            # Per-configuration jobs: these will wait/read the norm file and then run using best N₀.
            for N₀ in N₀_list
                for config_index in 1:nconfigs
                    randomseed = rand(Int64)

                    # pass N₀ = -1.0 as placeholder when we want the job to load the norm file
                    placeholder_N₀ = find_norm ? -1.0 : N₀

                    sbatch_command =
                        "sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index)
#SBATCH --output=/scratch/lappi/dana/slurm_out/runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index).out
#SBATCH --error=/scratch/lappi/dana/slurm_out/runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000

module load julia

# mode \"run\" will, if N₀ <= 0, wait for and read the normalization file written by the norm_only job.
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
"
                    write(io, sbatch_command)
                end
            end
        end
    end
end

println("Generated sbatch jobs in $output_file")