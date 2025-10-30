using Printf
using Serialization
using Random

if length(ARGS) < 1
    error("Usage: julia generatebatchjobsmaxnjobs.jl path/to/paramsfile.jl [njobsmax]")
end

params_file = ARGS[1]
njobs_max = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 300
if njobs_max <= 0
    error("njobsmax must be positive")
end

params_basename = split(basename(params_file), ".")[1]

include(params_file)
using .SimulationParams

outdir = joinpath("scripts", "sbatchjobs")
if !isdir(outdir)
    mkpath(outdir)
end

generated_files = String[]

for (set_index, params) in enumerate(SimulationParams.param_sets)
    m          = params.m
    nmax       = params.nmax
    sigma_list = params.sigma_list
    savepath   = params.savepath
    nconfigs   = params.nconfigs
    N₀_list    = params.N₀_list
    paramset   = params.paramset
    find_norm  = get(params, :find_norm, false)

    Ns = length(sigma_list)
    Nn0 = length(N₀_list)
    Nc = nconfigs

    total_jobs = Ns * Nn0 * Nc
    if total_jobs == 0
        @printf("Skipping param set %d (%s): no jobs\n", set_index, paramset)
        continue
    end

    # Write a single script containing all jobs for this param-set
    output_file = joinpath(outdir, "submitalljobs_$(params_basename)_set$(set_index).sh")
    open(output_file, "w") do io
        write(io, "#!/bin/bash\n\n")

        # Map linear job index -> (sigi, n0i, config_index)
        function job_indices_from_linear(j)
            t = j - 1
            cfg = (t % Nc) + 1
            t ÷= Nc
            n0i = (t % Nn0) + 1
            t ÷= Nn0
            sigi = t + 1
            return sigi, n0i, cfg
        end

        for j in 1:total_jobs
            sigi, n0i, config_index = job_indices_from_linear(j)
            sigma = sigma_list[sigi]
            N₀ = N₀_list[n0i]
            randomseed = rand(Int64)

            sbatch_command =
"""sbatch <<'EOF'
#!/bin/bash
#SBATCH --job-name=runshapefunctions_$(set_index)_$(paramset)_$(config_index)
#SBATCH --output=/scratch/lappi/dana/slurm_out/runshapefunctions_$(set_index)_$(paramset)_$(config_index).out
#SBATCH --error=/scratch/lappi/dana/slurm_out/runshapefunctions_$(set_index)_$(paramset)_$(config_index).err
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
    $(find_norm) \\
    $(N₀)
EOF

"""
            write(io, sbatch_command)
        end
    end
    push!(generated_files, output_file)
end

println("Generated $(length(generated_files)) sbatch job files:")
for f in generated_files
    println("  ", f)
end