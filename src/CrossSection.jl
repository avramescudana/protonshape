module CrossSection

"""
Packages
"""

using SpecialFunctions # Modified Bessel functions of second kind
using Symbolics # Symbolic calculation, partial derivatives
using MCIntegration # MC algorithms for high-dimensional integrals
using Distributions # Random numbers, Gaussian distributions
using Base.Threads # Multithreading

# using LinearAlgebra # Useful functions
# using SymbolicNumericIntegration # Symbolic integration
# using Nemo # Algebra package, requires for symbolic_solve
# using DelimitedFiles # Read data from text files
# using NLsolve # Numerically solve systems of equations
# using QuadGK # Numerically solve integrals
# using Statistics # Variance, standard deviation

"""
Parameters
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
    σ₀ = 29 * 2.56819, # [mb] 1 mb = 2.56819 GeV−2
    Λ = 0.28,
    x₀ = 4 * 10^(-5),
    Bₚ = 4, # [GeV^-2]
    N₀ = 50, # normalization 
    xₚ = 1.7 * 10^(-3), # corresponds to W = 75 GeV
    # xₚ = 9.6 * 10^(-4), # corresponds to W = 100 GeV
)

export params_gbw

# Conversion units
ħc = 0.197326 # [GeV*fm], convert GeV^-1 to fm
ħcinv = 5.068 # convert fm^-1 to GeV

params_mc = (
    rmin = 0.0,
    rmax = 5.0 * ħcinv, # rmax ≈ 1 fm, maximum dipole size < proton radius
    bmin = 0.0,
    bmax = 50 * ħcinv, # bmax ≈ 10 fm, maximum impact parameter
    zmin = 0.0,
    zmax = 1.0,
    θbmin = 0.0,
    θbmax = 2π,
    Δmin = 0.0,
    Δmax = 1.0,
    Δlen = 15,
    neval = 100000, # number of evaluations for MC integration
    niters = 10, # number of iterations for MC integration
)

export params_mc

diff_mode = "coh" # "coh" or "incoh"
dipole_mode = "GWB" # "GWB" or "CQ"

export diff_mode, dipole_mode

params_cq = (
    N₀ = 12, # normalization 
    Bqc = 3.3, # [GeV^-2]
    Bq = 0.7, # [GeV^-2]
    Nq = 3, # number of constituent quarks
    Nsamples = 50, # number of samples for bqc
)

export params_cq

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

include("dipole.jl") 
export Qₛ, T, gbwdipole
export sample_bqc, Tq, Tp, compute_Tp_grid

include("diffractive.jl")
export Agbw, Aqc, diffractive

end