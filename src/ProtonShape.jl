module ProtonShape

"""
Packages
"""

using SpecialFunctions # For Bessel 
using Symbolics # Symbolic calculation, partial derivatives
using MCIntegration # MC algorithms for high-dimensional integrals
using Distributions # Random numbers, Gaussian distributions
using Base.Threads # Multithreading
using Statistics # Variance, mean
using ProgressMeter # Progress bar
using Roots  # For root-finding
using Distributed # Multi-processing on single node
using ClusterManagers # Use by CSC for Slurm job submission
using JLD2 # Save to Julia file
using Serialization # Save to Julia file
using Interpolations 
using Optim

"""
Default parameters
"""

# Wavefunctions overlap parameters

αₑₘ = 1/137 # QED running coupling at Q²∼0y

params_wavefct = (
    Nₜ = 0.578, 
    Nₗ = 0.575,
    mf = 1.4, # [GeV] mass of charm quark, used in GBW fit
    R² = 2.3, # [GeV^-2]
    e = √(4π*αₑₘ), # QED coupling
    êf = 2, # effective charge, 2 for J/ψ
    Nc = 3, # number of colors
    Mᵥ = 3.097, # [GeV] mass of J/ψ
    δ = 1, # matched to other models, either 0 or 1
    Q² = 0, # data for J/ψ photoproduction at Q²=0
)

export params_wavefct

# GWB dipole model parameters 
# Fit to F₂ data including charm quarks 

params_gbw = (
    # σ₀ = 29 * 2.56819, # [mb] 1 mb = 2.56819 GeV−2
    # σ₀ = 29, # [mb] 
    Λ = 0.28,
    x₀ = 4 * 10^(-5),
    Bₚ = 4, # [GeV^-2]
    N₀ = 2, # normalization 
    xₚ = 1.7 * 10^(-3), # corresponds to W = 75 GeV
    # xₚ = 9.6 * 10^(-4), # corresponds to W = 100 GeV
)

export params_gbw

# Conversion units
ħc = 0.197326 # [GeV*fm], convert GeV^-1 to fm
ħcinv = 5.068 # convert fm^-1 to GeV

params_mc = (
    rmin = 0.0,
    rmax = 1.0 * ħcinv, # rmax ≈ 1 fm, maximum dipole size < proton radius
    bmin = 0.0,
    bmax = 10 * ħcinv, # bmax ≈ 10 fm, maximum impact parameter
    zmin = 0.0,
    zmax = 1.0,
    θbmin = 0.0,
    θbmax = 2π,
    Δmin = 0.0,
    Δmax = 1.0,
    Δlen = 20,
    # Δlen = 5,
    neval = 100000, # number of evaluations for MC integration
    niters = 10, # number of iterations for MC integration
    error_method = "jackknife", # error type for incoh, "standard", "halfsample" or "jackknife"
)

export params_mc

diff_mode = "coh" # "coh" or "incoh"
dipole_mode = "GWB" # "GWB" or "CQ"

export diff_mode, dipole_mode

params_cq = (
    N₀ = 1, # normalization 
    Bqc = 3.3, # [GeV^-2]
    Bq = 0.7, # [GeV^-2]
    Nq = 3, # number of constituent quarks
    Nsamples = 5, # number of samples for bqc
)

export params_cq

ħc = 0.197326 # [GeV*fm], convert GeV^-1 to fm

params_shape = (
    N₀ = 1.0, # normalization 
    α = 4.0, # gaussian radial function [GeV^-2]
    a = √8, # radius of the circular membrane [GeV^-1]
    σ = 5.0, # width of Gaussian distribution for amp, mean zero
    Nsamples = 1, # number of samples for amp
    # coeff_dict = Dict(), # dictionary with "(m,n) => amp" for the circular membrane
    type = "samem_multin", # type of sampling for the circular membrane
    # type = "samemn", # type of sampling for the circular membrane
    # current supported modes: samemn, samem_multin
    mn = (3,1), # (m,n) for the circular membrane
    nvals = 3, # number of n values for the circular membrane
    rotate = true, # random rotations in θb
    checktp = false, # check if Tp becomes negative
    replacetp = false, # replace negative Tp with zero
)

export params_shape

params_run = (
    run = "remote", # local for running on local machine, remote for cluster
    jobtype = "single", # job type for cluster, "array" or "single"
    # arrayindex = 1, # array index for cluster job
    savefile = true,
    run_threads = false,
    # savepath = "results/",
    savepath = "/scratch/lappi/dana/protonshape_out/mn_31/",
    outdir = "random/",
    # crosssec = "coh+incoh",
    # mode = "shapeamp",
    amp_dict = Dict{Tuple{Int,Int},Float64}(), # dictionary with "(m,n) => amp" for the circular membrane
)

export params_run

params_norm = (
    Δ₀ = 0.0, # fixed Δ value for normalization
    start = 1.0, # N₀ starting point for adaptive search
    step_factor = 1.6, # step factor for adaptive search
    min_N₀ = 0.001, # minimum N₀ for adaptive search
    max_N₀ = 100, # maximum N₀ for adaptive search
    max_expansions = 20, # maximum number of expansions for adaptive search
    unique_outdirs = false, # whether to use unique output directories for each run
    brent_reltol = 1e-3, # relative tolerance for Brent's method
    coherent_data_path = "/users/davrames/protonshape/data/jspi_coh_hera.txt", # path to coherent data for normalization
)

export params_norm

"""
Variables
"""

@variables r, b, θb, z, Δ
export r, b, θb, z, Δ

"""
Routines to extract coherent + incoherent dσ/dt
"""

include("wavefunction.jl") 
export ϕ, ΨᵥΨ

include("shape.jl")
export Tp_shape, shapedipole, sample_amp_dict_same_mn, sample_amp_dict_samem_multin, gaussenv, circmemb_2D, Tp_2D

include("dipole.jl") 
export Qₛ, T, gbwdipole
export sample_bqc, Tq, Tp, compute_Tp_grid

include("diffractive.jl")
export Agbw, Aqc, diffractive, compute_cross_sections

include("normalization.jl")
export chisq_for_N₀_at_Δ₀, find_best_N₀_at_Δ₀_adaptive

end