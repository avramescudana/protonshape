module SimulationParams

const param_sets = [
    (
        m = 3,
        nmax = 3,
        sigma_list = [5.0, 10.0, 15.0],
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_3_sigma_replacetp",  
        paramset = "m_3_nmax_3_sigma_replacetp", 
        nconfigs = 300,
        N₀_list = [1.0],
        find_norm = false,
        overrides = (
            params_shape = (replacetp = true, checktp = true),
        ),
    ),
]

export param_sets

end