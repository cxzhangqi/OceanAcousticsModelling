module OceanParameters

function bottom_loss_(θ₁)
	# Functions
	Wavelength(c, f) = c/f
	Wavenumber(λ) = 2π/λ
	Wavenumber(c, f) = Wavenumber(Wavelength(c, f))

	# Attenuation
	α(α_) = α_/(20log10(ℯ))
	ComplexSoundSpeed(cᵣ, α_, f) = cᵣ*(1 - im*α(α_)*cᵣ/2/π/f)

	# Signal
	f = 1e3

	# Water
	ρ₁ = 1e3
	c₁ = 1.5e3
	α₁ = 0
	λ₁ = Wavelength(c₁, f)
	k₁ = Wavenumber(c₁, f)
	ς₁ = ComplexSoundSpeed(c₁, α₁, f)

	# Sediment
	ρ₂ = 2e3
	cₚ = 1.6e3
	αₚ = 0.5
	cₛ = 0.
	αₛ = 0.
	λₚ = Wavelength(cₚ, f)
	kₚ = Wavenumber(cₚ, f)
	ςₚ = ComplexSoundSpeed(cₚ, αₚ, f)
	λₛ = Wavelength(cₛ, f)
	kₛ = Wavenumber(cₛ, f)
	ςₛ = ComplexSoundSpeed(cₛ, αₛ, f)

	# Bottom Loss
	kcosθ₁ = k₁*cos(θ₁)
	θₚ = invcos(kcosθ₁/kₚ)
	θₛ = invcos(kcosθ₁/kₛ)
	# 𝒵₁ = ρ₁*c₁/sin(θ₁)
	# 𝒵ₚ = ρ₂*cₚ/sin(θₚ)
	# 𝒵ₛ = ρ₂*cₛ/sin(θₛ)
	𝒵₁ = ρ₁*ς₁/sin(θ₁)
	𝒵ₚ = ρ₂*ςₚ/sin(θₚ)
	𝒵ₛ = ρ₂*ςₛ/sin(θₛ)
	𝒵_tot = 𝒵ₚ*cos(2θₛ)^2 + 𝒵ₛ*sin(2θₛ)^2
	ℛ = (𝒵_tot - 𝒵₁)/(𝒵_tot + 𝒵₁)
	BL = -10log10(abs(ℛ)^2)
end

"""

"""
function bottom_loss(θ₁; cₚ = 1.6e3, αₚ = 0.5, ρ₂ = 2e3, cₛ = 0.)
	# Functions
	Wavelength(c, f) = c/f
	Wavenumber(λ) = 2π/λ
	# Wavenumber(c, f) = Wavenumber(Wavelength(c, f))
	MyInv(val) = val == 0 ? Inf : 1/val
	Wavenumber(c, f) = 2π*f*MyInv(c)
	invcos(z) = π/2 + im*log(im*z + sqrt(1 - z^2))
	LossTangent(α_dbpwl) = α_dbpwl/(40π*log10(ℯ))
	ComplexSoundSpeed(cᵣ, δ) = cᵣ/(1 - im*δ)

	# Signal
	f = 1e3

	# Water
	ρ₁ = 1e3
	c₁ = 1.5e3
	α₁ = 0
	δ₁ = LossTangent(α₁)
	ς₁ = ComplexSoundSpeed(c₁, δ₁)
	λ₁ = Wavelength(ς₁, f)
	k₁ = Wavenumber(ς₁, f)

	# Sediment: Compressional
	δₚ = LossTangent(αₚ)
	ςₚ = ComplexSoundSpeed(cₚ, δₚ)
	λₚ = Wavelength(ςₚ, f)
	kₚ = Wavenumber(ςₚ, f)

	# Sediment: Shear
	αₛ = 0.
	δₛ = LossTangent(αₛ)
	ςₛ = ComplexSoundSpeed(cₛ, δₛ)
	λₛ = Wavelength(ςₛ, f)
	kₛ = Wavenumber(ςₛ, f)

	# Bottom Loss
	kcosθ₁ = k₁*cos(θ₁)
	θₚ = invcos(Complex(kcosθ₁/kₚ))
	θₛ = invcos(Complex(kcosθ₁/kₛ))
	𝒵₁ = ρ₁*ς₁/sin(θ₁)
	𝒵ₚ = ρ₂*ςₚ/sin(θₚ)
	𝒵ₛ = ρ₂*ςₛ/sin(θₛ)
	𝒵_tot = 𝒵ₚ*cos(2θₛ)^2 + 𝒵ₛ*sin(2θₛ)^2
	ℛ = (𝒵_tot - 𝒵₁)/(𝒵_tot + 𝒵₁)
	BL = -10log10(abs(ℛ)^2)

	return BL
end

end