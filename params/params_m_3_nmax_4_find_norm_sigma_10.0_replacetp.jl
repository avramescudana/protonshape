module SimulationParams

const param_sets = [
    (
        m = 3,
        nmax = 4,
        sigma_list = [10.0],
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_4_find_norm_sigma_10.0_replacetp",
        paramset = "m_3_nmax_4_find_norm_sigma_10.0_replacetpt", 
        nconfigs = 300,
        N₀_list = [0.6], # Placeholder, not used when find_norm=true
        find_norm = true,
        nsamples_norm = 20,  

        # Per-run overrides of ProtonShape params 
        overrides = (
            params_shape = (
                checktp   = true,
                replacetp = true,
            ),
            params_norm = (
                nsamples_norm = 20,   
            ),
        ),
    ),
]

export param_sets

end
