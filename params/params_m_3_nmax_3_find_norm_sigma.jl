module SimulationParams

const param_sets = [
    (
        m = 3,
        nmax = 3,
        sigma_list = [1.0],
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_3_find_norm_sigma",
        paramset = "m_3_nmax_3_find_norm_sigma", 
        nconfigs = 300,
        find_norm = true,
        N₀_list = [1.0], # Placeholder, not used when find_norm=true
    ),
]

export param_sets

end