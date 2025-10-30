module SimulationParams

const param_sets = [
    (
        m = 3,
        nmax = 3,
        sigma_list = [1.0, 5.0, 10.0],
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_3_find_norm_sigma",
        paramset = "m_3_nmax_3_find_norm_sigma", 
        nconfigs = 5,
        find_norm = true,
        N₀_list = [0.6], # Placeholder, not used when find_norm=true
        nsamples_norm = 1,
    ),
]

export param_sets

end