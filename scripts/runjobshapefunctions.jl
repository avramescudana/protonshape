include("../params/simulationparams.jl")
using .SimulationParams

function runshapefunction(config_index, randomseed)
    m = SimulationParams.m
    nmax = SimulationParams.nmax
    sigma_list = SimulationParams.sigma_list
    savepath = SimulationParams.savepath
    nconfigs = SimulationParams.nconfigs

    println("Running shape function with the following parameters:")
    println("Config Index: $config_index")
    println("Random Seed: $randomseed")
    println("m: $m")
    println("nmax: $nmax")
    println("sigma_list: $sigma_list")
    println("savepath: $savepath")
    println("nconfigs: $nconfigs")
end

for config_index in 1:SimulationParams.nconfigs
    randomseed = rand(Int)  
    runshapefunction(config_index, randomseed)
end