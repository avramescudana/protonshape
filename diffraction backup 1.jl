### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 89e17f34-ef6b-4537-b634-6f66f829489c
begin
	using SpecialFunctions # Modified Bessel functions of second kind
	using Symbolics # Symbolic calculation, partial derivatives
	using FFTW # Fast Fourier Transform, bindings to FFTW
	using LinearAlgebra # Useful LA tools
	using MCIntegration # MC algorithms for high-dimensional integrals
	# using SymbolicNumericIntegration # Symbolic integration
end

# ╔═╡ 9af71ee8-0317-11f0-102f-69a679def0dd
md"# Coherent and incoherent diffraction cross setion
Simple calculation of coferent and incoferent diffraction cross section for exclusive vector meson production. The calculation follows [arXiv:1607.01711](https://arxiv.org/abs/1607.01711). 

---
"

# ╔═╡ ec25d46c-ee6e-4b99-9a9c-3d36256223c0
md"### Import packages"

# ╔═╡ 47bd04c7-7db4-4de3-940f-3dfb7af1dd51
md"### Scattering amplitude
The scattering amplitude for diffractive vector meson production is given by

$$\begin{aligned}\mathcal{A}_{T,L}^{\gamma^*p\rightarrow Vp}(x_{\mathbb{P}}, Q^2, \boldsymbol{\Delta})=&\mathrm{i}\int \mathrm{d}^2 \boldsymbol{r}\int \mathrm{d}^2 \boldsymbol{b}\int\dfrac{\mathrm{d}z}{4\pi} (\Psi^*\Psi_V)_{T,L}(Q^2, \boldsymbol{r},z)\\ &\times \mathrm{e}^{-\mathrm{i}[\boldsymbol{b}-(1-z)\boldsymbol{r}]\cdot\boldsymbol{\Delta}}\dfrac{\mathrm{d}\sigma^p_{\mathrm{dip}}}{\mathrm{d}^2 \boldsymbol{b}}(\boldsymbol{b}, \boldsymbol{r}, x_{\mathbb{P}})\end{aligned}$$

The dipole-proton cross section $\sigma^p_\mathrm{dip}(\boldsymbol{b},\boldsymbol{r},x_\mathbb{P})$ is Fourier transformed to momentum space, with $\boldsymbol{\Delta}$ the Fourier conjugate to the center-of-mass of the dipole $\boldsymbol{b}-(1-z)\boldsymbol{r}$. Here $t=-\Delta^2$.

---
"

# ╔═╡ 893ba12a-c926-4cdc-912c-5c3a0f219dc2
md"#### Photon vector meson wavefunction overlap
The expressions are taken from [arXiv:hep-ph/0606272](https://arxiv.org/abs/hep-ph/0606272)
$$\begin{aligned}(\Psi^*_V\Psi)_T&=\hat{e}_f e\dfrac{N_c}{\pi z(1-z)}\left\{m_f^2 K_0(\epsilon r)\phi_T(r,z)-[z^2+(1-z)^2]\epsilon K_1(\epsilon r)\partial_r \phi_T(r,z)\right\}\\
(\Psi^*_V\Psi)_L&=\hat{e}_f e\dfrac{N_c}{\pi}2Qz(1-z)K_0(\epsilon r)\left[M_V\phi_L(r,z)+\delta \dfrac{m_f^2-\nabla_r^2}{M_vz(1-z)}\phi_L(r,z)\right]
\end{aligned}$$

---
"

# ╔═╡ 06a91492-5dc7-4fb5-814c-bfdd78beecc5
md"##### Boosted Gaussian model
$$\phi_{T,L}(r,z)=\mathcal{N}_{T,L}z(1-z)\exp\left(-\dfrac{m_f^2\mathcal{R}^2}{8z(1-z)}-\dfrac{2z(1-z)r^2}{\mathcal{R}^2}+\dfrac{m_f^2\mathcal{R}^2}{2}\right)$$

The parameters $\mathcal{N}_{T,L}$ and $\mathcal{R}$ are determined from the normalization and decay width conditions. 
"

# ╔═╡ 60264611-33c0-43ce-8e28-6687163cd669
md"##### Parameters"

# ╔═╡ 4cbad278-7963-47d8-8998-3e5d8a81637c
begin
	# wavefunction overlap parameters
	
	Nc = 3 # number of colors
	αₑₘ = 1/137 # QED running coupling at Q²∼0
	e = √(4π*αₑₘ) # QED coupling
	
	êf = 2 # effective charge, 2 for J/ψ
	
	# mf = 1.27 # [GeV] mass of charm quark
	mf = 1.5 # [GeV] mass of charm quark, used in GBW fit
	mₚ = 0.938 # [GeV] mass of proton
	
	δ = 1 # matched to other models, either 0 or 1
	Mᵥ = 3.1 # [GeV] mass of J/ψ

	# boosted gaussian wavefunction parameters
	# Table 2 from arXiv:hep-ph/0606272 for J/ψ

	Nₜ = 0.578 
	Nₗ = 0.575
	R² = 2.3 # [GeV^-2]

	Q² = 0 # data for J/ψ photoproduction at Q²=0
	W = 100 # [GeV]
end

# ╔═╡ 03f627f3-f72e-4f67-bc0c-822e491fd81d
function xp(t)
	nom = Mᵥ*Mᵥ + Q² - t
	den = W*W + Q² - mₚ*mₚ
	return nom/den
end

# ╔═╡ 96737eda-4f4b-427f-9644-a9b9540b360d
xp(t)

# ╔═╡ b0820c13-089b-485b-948f-103b20773c90
xₚ = xp(0)

# ╔═╡ d4f451d1-4a4b-4582-b4b9-9728bdebc25b
md"##### Functions"

# ╔═╡ 0b4680a1-b139-4720-87d4-62025e30feba
@variables r, z

# ╔═╡ 65f79988-11f5-4454-bac3-8cc680205019
ϵ(z, Q²) = sqrt(z * (1-z) * Q² + mf * mf)

# ╔═╡ 6b92005d-0098-42a0-9c7b-1b9fb36b97c5
function ϕ(r,z,λ)
	if λ=="T"
		norm = Nₜ
	elseif λ=="L"
		norm = Nₗ
	end
	factor = norm * z * (1-z)
	term1 = - mf * mf * R² / (8 * z * (1-z))
	term2 = - 2 * z * (1-z) * r * r / R²
	term3 = mf * mf * R² / 2
	term_exp = exp(term1 + term2 + term3)
	return factor * term_exp
end

# ╔═╡ a61eb470-b3f8-4739-8533-5476f8ff93bd
begin
	∂ᵣ = Differential(r)

	∂ᵣϕₜ(r,z) = expand_derivatives(∂ᵣ(ϕ(r,z,"T")))

	# ∇ᵣ² = 1/r * ∂ᵣ + ∂ᵣ * ∂ᵣ
	∇ᵣϕₗ(r,z) = 1/r * expand_derivatives(∂ᵣ(ϕ(r,z,"L"))) + expand_derivatives(∂ᵣ(∂ᵣ(ϕ(r,z,"L"))))
end

# ╔═╡ 69528f60-33c4-4c84-9e8d-c61c7cb4c599
function ΨᵥΨ(Q²,r,z,λ)
	factor = êf * e * Nc / π
	if λ=="T"
		factorλ = 1 / (z * (1-z))
		term1 = mf * mf * besselk(0, ϵ(z, Q²) * r) * ϕ(r,z,λ)
		term2 = - (z * z + (1-z) * (1-z)) * ϵ(z, Q²) * besselk(1, ϵ(z, Q²) * r) * ∂ᵣϕₜ(r,z)
	elseif λ=="L"
		factorλ = 2 * sqrt(Q²) * z * (1-z) * besselk(0, ϵ(z, Q²) * r)
		term1 = Mᵥ * ϕ(r,z,λ)
		term2 = δ / (Mᵥ * z * (1-z)) * (mf * mf * ϕ(r,z,λ) - ∇ᵣϕₗ(r,z))
	end
	return factor * factorλ * (term1 + term2) 
end

# ╔═╡ 927bce9d-b61e-4886-9775-404bf37e7c48
md"Symbolic expressions for $(\Psi_V\Psi^*)_{T,L}$"

# ╔═╡ 17bab355-00e5-4288-8955-2c5f5f66e77b
ΨᵥΨ(Q²,r,z,"T")

# ╔═╡ acef09fd-0934-41e4-bbda-fa23de9f6b33
ΨᵥΨ(Q²,r,z,"L")

# ╔═╡ 0096ee88-5586-4721-9792-51adf979fa71
md"
---

#### GBW dipole cross section

$$\sigma_{q\overline{q}}^\mathrm{GBW}(x,r)=\sigma_0 \left(1-\mathrm{e}^{-r^2 Q_s^2(x)/4}\right)$$

where $Q_s^2(x)=(x_0/x)^{\lambda_{\mathrm{GBW}}}$.

Let us introduce an impact parameter dependence 

$$\dfrac{\mathrm{d}\sigma_{q\overline{q}}}{\mathrm{d}^2\boldsymbol{b}}=\left(1-\mathrm{e}^{-\mathcal{N}r^2 Q_s^2(x,b)/4}\right)$$

where $Q_s^2(x,b)=(x_0/x)^{\lambda_{\mathrm{GBW}}}T(b)$ with a Gaussian thickness function

$$T(b)=\dfrac{1}{2\pi B_p}\mathrm{e}^{-\frac{b^2}{2B_p}}$$

Here the parameter $B_p$ is the proton width. Here $\mathcal{N}$ is a normalization factor chosen to fit $\mathrm{d}\sigma/\mathrm{d}t$.

---
"

# ╔═╡ 9d4c1960-73da-42d4-a64b-78b104387ad4
md"##### Parameters"

# ╔═╡ ad97b327-3cf1-4761-90cb-00a459be634a
begin
	# fit to F₂ data including charm quarks
	# mc = 1.5 GeV used in the fit
	
	σ₀ = 29 # [mb]
	Λ = 0.28
	x₀ = 4 * 10^(-5)

	Bₚ = 4 # [GeV^-2]
end

# ╔═╡ 9aec76c5-cfeb-4522-8fae-58842413be9e
md"##### Functions"

# ╔═╡ 8c104443-fbb5-4392-a714-e48ddce0768f
@variables t, b

# ╔═╡ 15e890db-7c65-4dba-abf9-298c42a7660e
# @variables xₚ

# ╔═╡ 43a5f98e-3d18-4e9a-9382-2b397bdcf3d5
# function xₚ(t)
# 	nom = Mᵥ*Mᵥ + Q² - t
# 	den = W*W + Q² - mₚ*mₚ
# 	return nom/den
# end

# ╔═╡ 34c2fd53-67ec-4fdc-b0fd-3490fceb3b0d
# xₚ(t)

# ╔═╡ 281fcffe-f9b2-4bad-8725-f6761544e112
function Qₛ(xₚ)
	return (x₀ / xₚ) ^ Λ
end

# ╔═╡ 23253186-7857-4051-a9e3-538cc6527419
function T(b)
	norm = 1 / (2π * Bₚ)
	term_exp = exp(- b * b / (2 * Bₚ * Bₚ))
	return norm * term_exp
end

# ╔═╡ cfe26cfa-cfa3-4a9a-aaf0-2a8ea0413670
T(b)

# ╔═╡ 54eb21fd-2a24-4b0d-b05b-78c19bdc729a
function σqq̅(xₚ,r)
	term_exp = 1 - exp(- r * r * Qₛ(xₚ) / 4)
	return σ₀ * term_exp
end

# ╔═╡ 837e8db4-effe-4131-bde3-972325a994ae
σqq̅(xₚ,r)

# ╔═╡ be723fcc-c84b-4757-aaaf-b07b5e418856
function dσqq̅db(xₚ,r,b)
	term_exp = 1 - exp(- r * r * Qₛ(xₚ) * T(b) / 4)
	return σ₀ * term_exp
end

# ╔═╡ e0447f20-2597-4236-a5bd-8ca1b1053213
dσqq̅db(xₚ,r,b)

# ╔═╡ 9c395345-b5c3-439c-a663-e815ecb287f6
md"
---

#### Coherent cross section
Coherent diffraction cross section

$$\dfrac{\mathrm{d}\sigma^{\gamma^*p\rightarrow Vp}_\mathrm{c}}{\mathrm{d}t}=\dfrac{1}{16\pi}\left|\langle\mathcal{A}_{T,L}^{\gamma^*p\rightarrow Vp}(x_{\mathbb{P}}, Q^2, \boldsymbol{\Delta})\rangle\right|^2$$

---
"

# ╔═╡ 26454bd8-f492-41a8-8ca1-820a0cc8fb5d
md"Extract 2D FFT (contains the integral over $\vec{b}$), perform the integral over $z$, then the integral over $\vec{r}$"

# ╔═╡ ef7db9b6-2a64-4d55-a0c0-21f29b9546da
md"#### Numerical grids"

# ╔═╡ 5e7c84f6-56b9-46f1-a1c6-5b9414a26d3e
begin
	ħc = 0.197326 # [GeV*fm], convert GeV^-1 to fm
	ħcinv = 5.068 # convert fm^-1 to GeV
	
	Nr, Nb, Nz = 100, 100, 50  # Number of points for r, b and z grids
	
	rmin, rmax = 0, 5 * ħcinv # rmax ≈ 1 fm, maximum dipole size < proton radius
	bmin, bmax = 0, 50 * ħcinv # bmax ≈ 10 fm, maximum impact parameter
	zmin, zmax = 0, 1 # z ∈ [0,1]

	θrmin, θrmax = 0, 2π # polar angle in r grid
	θzmin, θzmax = 0, 2π # polar angle in z grid
end

# ╔═╡ 223f2010-a5fb-4bdd-84b9-b2483a40a400
md"
---

#### Incoherent cross section
Incoherent diffraction cross section

$$\dfrac{\mathrm{d}\sigma^{\gamma^*p\rightarrow Vp}_\mathrm{inc}}{\mathrm{d}t}=\dfrac{1}{16\pi}\left(\left\langle\left|\mathcal{A}_{T,L}^{\gamma^*p\rightarrow Vp}(x_{\mathbb{P}}, Q^2, \boldsymbol{\Delta})\right|^2\right\rangle-\left|\langle\mathcal{A}_{T,L}^{\gamma^*p\rightarrow Vp}(x_{\mathbb{P}}, Q^2, \boldsymbol{\Delta})\rangle\right|^2\right)$$
"

# ╔═╡ Cell order:
# ╟─9af71ee8-0317-11f0-102f-69a679def0dd
# ╟─ec25d46c-ee6e-4b99-9a9c-3d36256223c0
# ╠═89e17f34-ef6b-4537-b634-6f66f829489c
# ╟─47bd04c7-7db4-4de3-940f-3dfb7af1dd51
# ╟─893ba12a-c926-4cdc-912c-5c3a0f219dc2
# ╟─06a91492-5dc7-4fb5-814c-bfdd78beecc5
# ╟─60264611-33c0-43ce-8e28-6687163cd669
# ╠═4cbad278-7963-47d8-8998-3e5d8a81637c
# ╠═03f627f3-f72e-4f67-bc0c-822e491fd81d
# ╠═96737eda-4f4b-427f-9644-a9b9540b360d
# ╠═b0820c13-089b-485b-948f-103b20773c90
# ╟─d4f451d1-4a4b-4582-b4b9-9728bdebc25b
# ╠═0b4680a1-b139-4720-87d4-62025e30feba
# ╠═a61eb470-b3f8-4739-8533-5476f8ff93bd
# ╠═65f79988-11f5-4454-bac3-8cc680205019
# ╠═6b92005d-0098-42a0-9c7b-1b9fb36b97c5
# ╠═69528f60-33c4-4c84-9e8d-c61c7cb4c599
# ╟─927bce9d-b61e-4886-9775-404bf37e7c48
# ╠═17bab355-00e5-4288-8955-2c5f5f66e77b
# ╠═acef09fd-0934-41e4-bbda-fa23de9f6b33
# ╟─0096ee88-5586-4721-9792-51adf979fa71
# ╟─9d4c1960-73da-42d4-a64b-78b104387ad4
# ╠═ad97b327-3cf1-4761-90cb-00a459be634a
# ╟─9aec76c5-cfeb-4522-8fae-58842413be9e
# ╠═8c104443-fbb5-4392-a714-e48ddce0768f
# ╠═15e890db-7c65-4dba-abf9-298c42a7660e
# ╠═43a5f98e-3d18-4e9a-9382-2b397bdcf3d5
# ╠═34c2fd53-67ec-4fdc-b0fd-3490fceb3b0d
# ╠═281fcffe-f9b2-4bad-8725-f6761544e112
# ╠═23253186-7857-4051-a9e3-538cc6527419
# ╠═cfe26cfa-cfa3-4a9a-aaf0-2a8ea0413670
# ╠═54eb21fd-2a24-4b0d-b05b-78c19bdc729a
# ╠═837e8db4-effe-4131-bde3-972325a994ae
# ╠═be723fcc-c84b-4757-aaaf-b07b5e418856
# ╠═e0447f20-2597-4236-a5bd-8ca1b1053213
# ╟─9c395345-b5c3-439c-a663-e815ecb287f6
# ╟─26454bd8-f492-41a8-8ca1-820a0cc8fb5d
# ╟─ef7db9b6-2a64-4d55-a0c0-21f29b9546da
# ╠═5e7c84f6-56b9-46f1-a1c6-5b9414a26d3e
# ╟─223f2010-a5fb-4bdd-84b9-b2483a40a400
