module SimulationParams

const param_sets = [
    (
        m = 3,
        nmax = 3,
        sigma_list = [1.0],
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_3_find_norm_sigma_replacetp_test",
        paramset = "m_3_nmax_3_find_norm_sigma_replacetp_test", 
        nconfigs = 3,
        N₀_list = [0.6], # Placeholder, not used when find_norm=true
        find_norm = true,
        nsamples_norm = 2,  

        # Per-run overrides of ProtonShape params 
        overrides = (
            params_shape = (
                checktp   = true,
                replacetp = true,
            ),
            # params_norm = (
            #     nsamples_norm = 2,   
            #     ngrid = 3,
            # ),
            # params_mc = (
            #     Δlen = 2,
            # ),
        ),
    ),
]

export param_sets

end
