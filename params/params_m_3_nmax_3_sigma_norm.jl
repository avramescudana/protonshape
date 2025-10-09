module SimulationParams

const param_sets = [
    (
        m = 3,
        nmax = 3,
        sigma_list = [0, 5, 10, 15],
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_3_fine_sigma",  
        paramset = "m_3_nmax_3_sigma_norm", 
        nconfigs = 300,
        N₀_list = [0.1, 1.0, 10.0],
    ),
]

export param_sets

end