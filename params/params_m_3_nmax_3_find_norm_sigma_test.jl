module SimulationParams

const param_sets = [
    (
        m = 3,
        nmax = 3,
        sigma_list = [1.0],
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_3_find_norm_sigma_test",
        paramset = "m_3_nmax_3_find_norm_sigma_test", 
        nconfigs = 3,
        N₀_list = [0.6], # Placeholder, not used when find_norm=true
        find_norm = true,
        nsamples_norm = 2,
    ),
]

export param_sets

end