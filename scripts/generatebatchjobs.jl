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
        paramset = params.paramset
        find_norm  = params.find_norm

        sigma_values_str = join(sigma_list, " ")

        for sigma in sigma_list
            for N₀ in N₀_list
                for config_index in 1:nconfigs
                    randomseed = rand(Int64)

                    sbatch_command =
                        "sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=runshapefunctions_$(set_index)_$(paramset)_$config_index
#SBATCH --output=/scratch/lappi/dana/slurm_out/runshapefunctions_$(set_index)_$(paramset)_$config_index.out
#SBATCH --error=/scratch/lappi/dana/slurm_out/runshapefunctions_$(set_index)_$(paramset)_$config_index.err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000

module load julia

julia --project=. scripts/runshapefunctions.jl \\
    $config_index \\
    $nconfigs \\
    $randomseed \\
    $m \\
    $nmax \\
    $savepath \\
    $sigma \\
    $find_norm \\
    $N₀ 
EOF
"
                    write(io, sbatch_command)
                end
            end
        end
    end
end

println("Generated sbatch jobs in $output_file")
