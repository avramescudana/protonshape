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

# determine full path to sbatch once
sbatch_cmd = chomp(read(`which sbatch`, String))

open(output_file, "w") do io
    write(io, "#!/bin/bash\n\n")
    write(io, "set -euo pipefail\n\n")

    for (set_index, params) in enumerate(SimulationParams.param_sets)
        m             = params.m
        nmax          = params.nmax
        sigma_list    = params.sigma_list
        savepath      = params.savepath
        nconfigs      = params.nconfigs
        N₀_list       = params.N₀_list
        paramset      = params.paramset
        find_norm     = params.find_norm
        nsamples_norm = get(params, :nsamples_norm, 20)  # <- set to 4 in your param file

        # ensure savepath + slurm_out
        mkpath(savepath)
        slurm_out_dir = joinpath(savepath, "slurm_out")
        mkpath(slurm_out_dir)

        overrides = get(params, :overrides, nothing)
        if overrides !== nothing
            override_file = joinpath(savepath, "params_override.jl2")
            open(override_file, "w") do f
                write(f, "const params_override = $(repr(overrides))\n")
            end
        end

        for sigma in sigma_list
            sigma_sanit    = replace(string(sigma), "." => "_")
            paramset_sanit = replace(string(paramset), r"[^A-Za-z0-9_]" => "_")

            if mode == "csc"
                if find_norm
                    norm_randomseed = rand(Int64)

                    # 1) grid job
                    sbatch_grid_command = """
jobid_grid_$(set_index)_$(paramset_sanit)_$(sigma_sanit)=\$(
""" * sbatch_cmd * """ --parsable <<'GRID_EOF'
#!/bin/bash
#SBATCH --job-name=grid_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=$(slurm_out_dir)/grid_$(set_index)_$(paramset)_$(sigma).out
#SBATCH --error=$(slurm_out_dir)/grid_$(set_index)_$(paramset)_$(sigma).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1000

module load julia

julia --project=. scripts/runshapefunctions.jl \\
    1 \\
    $(nsamples_norm) \\
    $(norm_randomseed) \\
    $(m) \\
    $(nmax) \\
    $(savepath) \\
    $(sigma) \\
    true \\
    -1.0 \\
    $(paramset) \\
    norm_grid
GRID_EOF
)
"""
                    write(io, sbatch_grid_command * "\n")
                end

                # snippet to wait on norm_done marker in run jobs
                if find_norm
                    wait_snippet = "norm_done=\"" * savepath * "/norm/" * string(paramset) * "_sigma_" * string(sigma) * "_norm_done\"\n" *
                                   "while [ ! -f \"\$norm_done\" ]; do sleep 10; done\n"
                else
                    wait_snippet = ""
                end

                if find_norm
                    # 2) submitter job
                    sbatch_submitter_command = """
# resolve the dynamic grid job id and validate it before using as dependency
varname=jobid_grid_$(set_index)_$(paramset_sanit)_$(sigma_sanit)
gridid=\${!varname}

gridid=\$(echo "\$gridid" | tr -d '[:space:]')
gridid=\${gridid%%;*}

if [ -z "\$gridid" ]; then
  echo "submitter: jobid_grid not set (varname=\$varname)" >&2
  exit 1
fi
echo "submitter: using grid job id = \$gridid"

jobid_submitter_$(set_index)_$(paramset_sanit)_$(sigma_sanit)=\$(
""" * sbatch_cmd * """ --parsable --dependency=afterok:\$gridid <<'SUBMITTER_EOF'
#!/bin/bash
#SBATCH --job-name=submit_array_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=$(slurm_out_dir)/submit_array_$(set_index)_$(paramset)_$(sigma).out
#SBATCH --error=$(slurm_out_dir)/submit_array_$(set_index)_$(paramset)_$(sigma).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1000

set -euo pipefail

grid_txt="$(savepath)/norm/$(paramset)_sigma_$(sigma)_N0grid.txt"
if [ ! -f "\$grid_txt" ]; then
  echo "submitter: grid file not found: \$grid_txt" >&2
  exit 1
fi

N=\$(wc -w < "\$grid_txt" | tr -d ' ')
if [ -z "\$N" ] || [ "\$N" -lt 1 ]; then
  echo "submitter: invalid N from grid: \$N" >&2
  exit 1
fi

jobid_array=\$(
""" * sbatch_cmd * """ --parsable --array=1-\${N} <<'ARRAY_EOF'
#!/bin/bash
#SBATCH --job-name=norm_point_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=$(slurm_out_dir)/norm_point_$(set_index)_$(paramset)_$(sigma).%A_%a.out
#SBATCH --error=$(slurm_out_dir)/norm_point_$(set_index)_$(paramset)_$(sigma).%A_%a.err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4000

module load julia

grid_txt="$(savepath)/norm/$(paramset)_sigma_$(sigma)_N0grid.txt"
N0vals=(\$(cat "\$grid_txt"))
N0=\${N0vals[\$((\$SLURM_ARRAY_TASK_ID-1))]}

# arrayindex = 1 for normalization; N0 is passed separately
julia --project=. scripts/runshapefunctions.jl \\
    1 \\
    $(nsamples_norm) \\
    $(norm_randomseed) \\
    $(m) \\
    $(nmax) \\
    $(savepath) \\
    $(sigma) \\
    true \\
    \${N0} \\
    $(paramset) \\
    norm_point
ARRAY_EOF
)

jobid_array=\$(echo "\$jobid_array" | tr -d '[:space:]')
jobid_array=\${jobid_array%%;*}
if [ -z "\$jobid_array" ]; then
  echo "submitter: failed to submit array (empty jobid)" >&2
  exit 1
fi

jobid_collect=\$(
""" * sbatch_cmd * """ --parsable --dependency=afterok:\${jobid_array} <<'COLLECT_EOF'
#!/bin/bash
#SBATCH --job-name=norm_collect_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=$(slurm_out_dir)/norm_collect_$(set_index)_$(paramset)_$(sigma).out
#SBATCH --error=$(slurm_out_dir)/norm_collect_$(set_index)_$(paramset)_$(sigma).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=2000

module load julia

grid_txt="$(savepath)/norm/$(paramset)_sigma_$(sigma)_N0grid.txt"
julia --project=. scripts/runshapefunctions.jl \\
    1 \\
    $(nsamples_norm) \\
    $(norm_randomseed) \\
    $(m) \\
    $(nmax) \\
    $(savepath) \\
    $(sigma) \\
    true \\
    -1.0 \\
    $(paramset) \\
    norm_collect \\
    "\$grid_txt"

mkdir -p "$(savepath)/norm"
touch "$(savepath)/norm/$(paramset)_sigma_$(sigma)_norm_done"
COLLECT_EOF
)

echo "submitter: submitted array \$jobid_array and collector \$jobid_collect"
SUBMITTER_EOF
)

jobid_submitter_$(set_index)_$(paramset_sanit)_$(sigma_sanit)=\$(echo "\$jobid_submitter_$(set_index)_$(paramset_sanit)_$(sigma_sanit)" | tr -d '[:space:]')
jobid_submitter_$(set_index)_$(paramset_sanit)_$(sigma_sanit)=\${jobid_submitter_$(set_index)_$(paramset_sanit)_$(sigma_sanit)%%;*}
if [ -z "\$jobid_submitter_$(set_index)_$(paramset_sanit)_$(sigma_sanit)" ]; then
  echo "failed to capture submitter job id for set=$(set_index) paramset=$(paramset) sigma=$(sigma)" >&2
  exit 1
fi
echo "top-level: submitter job id = \$jobid_submitter_$(set_index)_$(paramset_sanit)_$(sigma_sanit)"
"""
                    write(io, sbatch_submitter_command * "\n")
                end

                # 3) run jobs (physics) depend on submitter and wait for norm_done marker
                for N₀ in N₀_list
                    for config_index in 1:nconfigs
                        randomseed     = rand(Int64)
                        placeholder_N₀ = find_norm ? -1.0 : N₀

                        dependency_arg = ""
                        if find_norm
                            jobid_varname  = "jobid_submitter_$(set_index)_$(paramset_sanit)_$(sigma_sanit)"
                            dependency_arg = "--dependency=afterok:\${" * jobid_varname * "}"
                        end

                        dep_prefix = isempty(dependency_arg) ? "" : dependency_arg * " "

                        sbatch_command = """
""" * sbatch_cmd * """ """ * dep_prefix * """<<'EOF'
#!/bin/bash
#SBATCH --job-name=runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index)
#SBATCH --output=$(slurm_out_dir)/runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index).out
#SBATCH --error=$(slurm_out_dir)/runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000

module load julia

$(wait_snippet)julia --project=. scripts/runshapefunctions.jl \\
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
                # local mode
                if find_norm
                    norm_randomseed = rand(Int64)
                    cmd = @sprintf("echo 'Running norm_only: set=%d sigma=%s'\nmkdir -p %s\njulia --project=. scripts/runshapefunctions.jl %d %d %d %d %d %s %g true %g %s norm_only\n\n",
                        set_index, string(sigma), savepath,
                        1, nsamples_norm, norm_randomseed, m, nmax, savepath, sigma, -1.0, paramset)
                    write(io, cmd)
                end

                for N₀ in N₀_list
                    for config_index in 1:nconfigs
                        randomseed     = rand(Int64)
                        placeholder_N₀ = find_norm ? -1.0 : N₀
                        cmd = @sprintf("echo 'Running job: set=%d paramset=%s sigma=%g config=%d'\nmkdir -p %s\njulia --project=. scripts/runshapefunctions.jl %d %d %d %d %d %s %g false %g %s run\n\n",
                            set_index, paramset, sigma, config_index, savepath,
                            config_index, nconfigs, randomseed, m, nmax, savepath, sigma, placeholder_N₀, paramset)
                        write(io, cmd)
                    end
                end
            end
        end
    end
end

println("Generated $mode jobs script: $output_file")
