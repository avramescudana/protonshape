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

    # Determine grouping so we produce at most `njobs_max` sbatch files (parts).
    parts_count = min(njobs_max, total_jobs)
    group_size = ceil(Int, total_jobs / parts_count)

    # Build linear list of all (sigma_index, N0_index, config_index) in a deterministic order
    jobs = Vector{NamedTuple{(:sigi, :n0i, :config_index, :sigma, :N0, :randomseed),Tuple{Int,Int,Int,Float64,Float64,Int64}}}()
    for sigi in 1:Ns, n0i in 1:Nn0, cfg in 1:Nc
        push!(jobs, (sigi=sigi, n0i=n0i, config_index=cfg,
                     sigma=sigma_list[sigi], N0=N₀_list[n0i], randomseed=rand(Int64)))
    end

    # Partition into parts and write one sbatch batch script per part.
    nparts = ceil(Int, length(jobs) / group_size)
    for part_idx in 1:nparts
        start_i = (part_idx - 1) * group_size + 1
        end_i = min(length(jobs), part_idx * group_size)
        part_jobs = jobs[start_i:end_i]

        output_file = joinpath(outdir, "submitalljobs_$(params_basename)_set$(set_index)_part$(part_idx).sh")
        open(output_file, "w") do io
            write(io, "#!/bin/bash\n")
            write(io, "#SBATCH --job-name=runshapefunctions_$(set_index)_$(paramset)_part$(part_idx)\n")
            write(io, "#SBATCH --output=/scratch/lappi/dana/slurm_out/runshapefunctions_$(set_index)_$(paramset)_part$(part_idx).out\n")
            write(io, "#SBATCH --error=/scratch/lappi/dana/slurm_out/runshapefunctions_$(set_index)_$(paramset)_part$(part_idx).err\n")
            write(io, "#SBATCH --account=lappi\n")
            write(io, "#SBATCH --partition=small\n")
            write(io, "#SBATCH --time=24:00:00\n")
            write(io, "#SBATCH --mem-per-cpu=4000\n\n")

            write(io, "module load julia\n\n")
            write(io, "set -euo pipefail\n\n")

            # If this param-set requested a single normalization job per sigma, create norm_only runs
            # for unique sigmas that appear in this part. Place them at the start of the script.
            if find_norm
                sigmas_in_part = unique(p.sigma for p in part_jobs)
                for σ in sigmas_in_part
                    norm_randomseed = rand(Int64)
                    write(io, "# Find normalization for sigma=$(σ)\n")
                    # call with nsamples_norm handled by params in the params file (runshapefunctions will read savepath)
                    write(io, @sprintf("julia --project=. scripts/runshapefunctions.jl %d %d %d %d %d \"%s\" %0.16g true -1.0 \"%s\" norm_only\n",
                                       1,  # placeholder arrayindex for norm-only call
                                       1,  # nconfigs (not used there)
                                       norm_randomseed,
                                       m,
                                       nmax,
                                       savepath,
                                       σ,
                                       paramset))
                    write(io, "\n")
                end
            end

            # Sequentially run all configurations in this part
            for pj in part_jobs
                N0_arg = find_norm ? -1.0 : pj.N0
                write(io, @sprintf("# config: sigma=%.6g N0=%.6g config_index=%d\n", pj.sigma, pj.N0, pj.config_index))
                write(io, @sprintf("julia --project=. scripts/runshapefunctions.jl %d %d %d %d %d \"%s\" %0.16g %s %0.16g \"%s\" run\n",
                                   pj.config_index,
                                   Nc,
                                   pj.randomseed,
                                   m,
                                   nmax,
                                   savepath,
                                   pj.sigma,
                                   string(find_norm),
                                   N0_arg,
                                   paramset))
                write(io, "\n")
            end
        end
        # make file executable (so user can sbatch it)
        try
            run(`chmod +x $output_file`)
        catch
            # ignore chmod errors on platforms without `chmod` available
        end
        push!(generated_files, output_file)
    end
end

# Create a single master submit script that sbatches each generated part/script.
# This is the one you run to submit all parts (it contains sbatch calls for every generated part).
if !isempty(generated_files)
    master_file = joinpath(outdir, "submit_all_$(params_basename).sh")
    open(master_file, "w") do io
        write(io, "#!/bin/bash\nset -euo pipefail\n\n")
        write(io, "# Master submit script generated by generatebatchjobsmaxnjobs.jl\n\n")
        for f in generated_files
            write(io, "echo 'Submitting: $(basename(f))'\n")
            write(io, @sprintf("sbatch \"%s\"\n\n", f))
        end
    end
    try
        run(`chmod +x $master_file`)
    catch
        # ignore chmod failures
    end
    println("Generated $(length(generated_files)) sbatch job files:")
    for f in generated_files
        println("  ", f)
    end
    println("Master submit script:")
    println("  ", master_file)
else
    println("No sbatch job files generated.")
end