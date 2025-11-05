using Printf
using Serialization
using Random

if length(ARGS) < 1
    error("Usage: julia generatebatchjobs.jl path/to/paramsfile.jl [mode]\n  mode: csc (default) | local")
end

params_file = ARGS[1]
mode        = length(ARGS) >= 2 ? ARGS[2] : "csc"
mode ∉ ("csc", "local") && error("Unknown mode: $mode — must be 'csc' or 'local'")

params_basename = split(basename(params_file), ".")[1]
output_dir      = joinpath("scripts", "sbatchjobs")
mkpath(output_dir)

output_file =
    mode == "csc" ?
        joinpath(output_dir, "submitalljobs_$(params_basename).sh") :
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
        nsamples_norm = get(params, :nsamples_norm, 20)

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

                    # 1) GRID JOB: compute N0 grid (norm_grid)
                    sbatch_grid_command = """
# === GRID job for set=$(set_index), paramset=$(paramset), sigma=$(sigma) ===
jobid_grid_$(set_index)_$(paramset_sanit)_$(sigma_sanit)=\$(
$sbatch_cmd --parsable <<'GRID_EOF'
#!/bin/bash
#SBATCH --job-name=grid_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=$(slurm_out_dir)/grid_$(set_index)_$(paramset)_$(sigma).out
#SBATCH --error=$(slurm_out_dir)/grid_$(set_index)_$(paramset)_$(sigma).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1000

module load julia

echo "arrayindex = 1"
echo "savepath   = $(savepath)"
echo "sigma      = $(sigma)"

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
echo "Submitted GRID job: \$jobid_grid_$(set_index)_$(paramset_sanit)_$(sigma_sanit)"

"""
                    write(io, sbatch_grid_command * "\n")

                    # 2) SUBMITTER JOB:
                    #    - submits norm_point array (one task per N0)
                    #    - submits norm_collect (depend on norm_point array)
                    #    - submits physics runshape array (depend on norm_collect)
                    sbatch_submitter_command = """
# === SUBMITTER job for set=$(set_index), paramset=$(paramset), sigma=$(sigma) ===
jobid_submitter_$(set_index)_$(paramset_sanit)_$(sigma_sanit)=\$(
$sbatch_cmd --parsable --dependency=afterok:\$jobid_grid_$(set_index)_$(paramset_sanit)_$(sigma_sanit) <<'SUBMITTER_EOF'
#!/bin/bash
#SBATCH --job-name=submit_norm_and_run_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=$(slurm_out_dir)/submit_array_$(set_index)_$(paramset)_$(sigma).out
#SBATCH --error=$(slurm_out_dir)/submit_array_$(set_index)_$(paramset)_$(sigma).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1000

set -euo pipefail

echo "SUBMITTER: starting for set=$(set_index) paramset=$(paramset) sigma=$(sigma)"

norm_dir="$(savepath)/norm"
mkdir -p "\$norm_dir"

grid_txt="\${norm_dir}/$(paramset)_sigma_$(sigma)_N0grid.txt"
if [ ! -f "\$grid_txt" ]; then
  echo "SUBMITTER: grid file not found: \$grid_txt" >&2
  exit 1
fi

N=\$(wc -w < "\$grid_txt" | tr -d ' ')
if [ -z "\$N" ] || [ "\$N" -lt 1 ]; then
  echo "SUBMITTER: invalid N from grid: \$N" >&2
  exit 1
fi
echo "SUBMITTER: N0 grid has \$N points"

# --- (a) norm_point array: one task per N0 grid point ---
jobid_array=\$(
$sbatch_cmd --parsable --array=1-\${N} <<'NORM_ARRAY_EOF'
#!/bin/bash
#SBATCH --job-name=norm_point_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=$(slurm_out_dir)/norm_point_$(set_index)_$(paramset)_$(sigma).%A_%a.out
#SBATCH --error=$(slurm_out_dir)/norm_point_$(set_index)_$(paramset)_$(sigma).%A_%a.err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=4000

module load julia

grid_txt="$(savepath)/norm/$(paramset)_sigma_$(sigma)_N0grid.txt"
N0vals=(\$(cat "\$grid_txt"))
N0=\${N0vals[\$((SLURM_ARRAY_TASK_ID-1))]}

echo "norm_point task: SLURM_ARRAY_TASK_ID=\$SLURM_ARRAY_TASK_ID, N0=\$N0"

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
    \$N0 \\
    $(paramset) \\
    norm_point
NORM_ARRAY_EOF
)

jobid_array=\$(echo "\$jobid_array" | tr -d '[:space:]')
jobid_array=\${jobid_array%%;*}
if [ -z "\$jobid_array" ]; then
  echo "SUBMITTER: failed to submit norm_point array (empty jobid)" >&2
  exit 1
fi
echo "SUBMITTER: norm_point array jobid = \$jobid_array"

# --- (b) norm_collect: depends on array ---
jobid_collect=\$(
$sbatch_cmd --parsable --dependency=afterok:\${jobid_array} <<'COLLECT_EOF'
#!/bin/bash
#SBATCH --job-name=norm_collect_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=$(slurm_out_dir)/norm_collect_$(set_index)_$(paramset)_$(sigma).out
#SBATCH --error=$(slurm_out_dir)/norm_collect_$(set_index)_$(paramset)_$(sigma).err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=01:00:00
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

COLLECT_EOF
)

jobid_collect=\$(echo "\$jobid_collect" | tr -d '[:space:]')
jobid_collect=\${jobid_collect%%;*}
if [ -z "\$jobid_collect" ]; then
  echo "SUBMITTER: failed to submit norm_collect (empty jobid)" >&2
  exit 1
fi
echo "SUBMITTER: norm_collect jobid = \$jobid_collect"

# --- (c) physics run array: depends on norm_collect ---
jobid_run=\$(
$sbatch_cmd --parsable --dependency=afterok:\${jobid_collect} --array=1-$(nconfigs) <<'RUN_ARRAY_EOF'
#!/bin/bash
#SBATCH --job-name=runshape_$(set_index)_$(paramset)_$(sigma)
#SBATCH --output=$(slurm_out_dir)/runshapefunctions_$(set_index)_$(paramset)_$(sigma).%A_%a.out
#SBATCH --error=$(slurm_out_dir)/runshapefunctions_$(set_index)_$(paramset)_$(sigma).%A_%a.err
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000

module load julia

cfg=\$SLURM_ARRAY_TASK_ID
randomseed=\$((123456 + cfg))

echo "runshape array task: cfg=\$cfg randomseed=\$randomseed"

julia --project=. scripts/runshapefunctions.jl \\
    \$cfg \\
    $(nconfigs) \\
    \$randomseed \\
    $(m) \\
    $(nmax) \\
    $(savepath) \\
    $(sigma) \\
    false \\
    -1.0 \\
    $(paramset) \\
    run
RUN_ARRAY_EOF
)

jobid_run=\$(echo "\$jobid_run" | tr -d '[:space:]')
jobid_run=\${jobid_run%%;*}
if [ -z "\$jobid_run" ]; then
  echo "SUBMITTER: failed to submit runshape array (empty jobid)" >&2
  exit 1
fi
echo "SUBMITTER: runshape array jobid = \$jobid_run"

echo "SUBMITTER: done. norm_point=\$jobid_array norm_collect=\$jobid_collect runshape=\$jobid_run"

SUBMITTER_EOF
)
echo "Submitted SUBMITTER job: \$jobid_submitter_$(set_index)_$(paramset_sanit)_$(sigma_sanit)"

"""
                    write(io, sbatch_submitter_command * "\n")

                else
                    # ---- CSC mode WITHOUT normalization (no dependencies) ----
                    for N₀ in N₀_list
                        for config_index in 1:nconfigs
                            randomseed     = rand(Int64)
                            sbatch_command = """
$sbatch_cmd <<'EOF'
#!/bin/bash
#SBATCH --job-name=runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index)
#SBATCH --output=$(slurm_out_dir)/runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index).out
#SBATCH --error=$(slurm_out_dir)/runshapefunctions_$(set_index)_$(paramset)_$(sigma)_$(config_index).err
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
    $(N₀) \\
    $(paramset) \\
    run
EOF

"""
                            write(io, sbatch_command * "\n")
                        end
                    end
                end
            else
                # -----------------------------
                # LOCAL mode (no Slurm)
                # -----------------------------
                if find_norm
                    norm_randomseed = rand(Int64)
                    cmd = @sprintf("""
echo 'Running norm_only: set=%d sigma=%s'
mkdir -p %s
julia --project=. scripts/runshapefunctions.jl %d %d %d %d %d %s %g true %g %s norm_only

""",
                        set_index, string(sigma), savepath,
                        1, nsamples_norm, norm_randomseed, m, nmax, savepath, sigma, -1.0, paramset)
                    write(io, cmd)
                end

                for N₀ in N₀_list
                    for config_index in 1:nconfigs
                        randomseed     = rand(Int64)
                        placeholder_N₀ = find_norm ? -1.0 : N₀
                        cmd = @sprintf("""
echo 'Running job: set=%d paramset=%s sigma=%g config=%d'
mkdir -p %s
julia --project=. scripts/runshapefunctions.jl %d %d %d %d %d %s %g false %g %s run

""",
                            set_index, paramset, sigma, config_index, savepath,
                            config_index, nconfigs, randomseed, m, nmax, savepath, sigma, placeholder_N₀, paramset)
                        write(io, cmd)
                    end
                end
            end
        end
    end
end

# make generated script executable
run(`chmod +x $output_file`)

println("Generated $mode jobs script: $output_file")
