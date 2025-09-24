using Printf
using Serialization

include("../params/simulationparams.jl")

open("submitalljobs.sh", "w") do io
    write(io, "#!/bin/bash\n\n")

    for config_index in 1:nconfigs
        randomseed = rand(Int64)
        sigma_values_str = join(sigma_list, ",")  

        sbatch_command = """
        sbatch <<EOF
        #!/bin/bash
        #SBATCH --job-name=runshapefunctions_$config_index
        #SBATCH --output=/scratch/lappi/dana/slurm_out/runshapefunctions_$config_index.out
        #SBATCH --error=/scratch/lappi/dana/slurm_out/runshapefunctions_$config_index.err
        #SBATCH --account=lappi
        #SBATCH --partition=small
        #SBATCH --time=24:00:00
        #SBATCH --mem-per-cpu=4000

        module load julia

        julia --project=. scripts/runshapefunctions.jl \
            "$config_index" \
            "$nconfigs" \
            "$randomseed" \
            "$m" \
            "$nmax" \
            "$savepath" \
            [$sigma_values_str]
        EOF
        """

        write(io, sbatch_command)
    end
end

println("Generated sbatch jobs in submitalljobs.sh")