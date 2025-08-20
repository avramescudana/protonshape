module SimulationParams

const param_sets = [
    (
        m = 4,
        nmax = 2,
        sigma_list = [8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10, 10.25, 10.5, 10.75, 11],
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_2_fine_sigma",  
        nconfigs = 300,
    ),
    (
        m = 4,
        nmax = 3,
        sigma_list = [8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10, 10.25, 10.5, 10.75, 11], 
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_3_fine_sigma",  
        nconfigs = 300,
    ),
    (
        m = 4,
        nmax = 4,
        sigma_list = [8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10, 10.25, 10.5, 10.75, 11], 
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_4_fine_sigma",  
        nconfigs = 300,
    ),
    (
        m = 4,
        nmax = 5,
        sigma_list = [8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10, 10.25, 10.5, 10.75, 11], 
        savepath = "/scratch/lappi/dana/protonshape_out/m_3_nmax_5_fine_sigma",  
        nconfigs = 300,
    ),
]

export param_sets

end
