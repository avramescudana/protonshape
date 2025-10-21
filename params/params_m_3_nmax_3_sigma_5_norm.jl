module SimulationParams

const param_sets = [
    (
        m = 3,
        nmax = 3,
        sigma_list = [5.0],
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_3_sigma_5_norm",  
        paramset = "m_3_nmax_3_sigma_5_norm", 
        nconfigs = 100,
        N₀_list = [0.1, 1.0, 10.0],
    ),
]

export param_sets

end