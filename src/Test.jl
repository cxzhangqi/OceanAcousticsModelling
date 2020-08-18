using Plots

include("OceanParameters.jl")

θ₁ = π/4
@show BL = OceanParameters.bottom_loss(θ₁)

θ₁ = deg2rad.(range(0, 90, length = 100))

cₚ = [1550 1600 1800]
Numθ₁ = length(θ₁)
Numcₚ = length(cₚ)
BL = zeros(AbstractFloat, (Numθ₁, Numcₚ))
for nθ₁ = 1:length(θ₁)
	for ncₚ = 1:length(cₚ)
		BL[nθ₁, ncₚ] = OceanParameters.bottom_loss(θ₁[nθ₁]; cₚ = cₚ[ncₚ])
	end
end

plot(rad2deg.(θ₁), BL)

αₚ = [1. 0.5 0]
Numθ₁ = length(θ₁)
Numαₚ = length(αₚ)
BL = zeros(AbstractFloat, (Numθ₁, Numαₚ))
for nθ₁ = 1:length(θ₁)
	for nαₚ = 1:length(αₚ)
		BL[nθ₁, nαₚ] = OceanParameters.bottom_loss(θ₁[nθ₁]; αₚ = αₚ[nαₚ])
	end
end

plot(rad2deg.(θ₁), BL)

ρ₂ = [1.5e3 2e3 2.5e3]
Numρ₂ = length(ρ₂)
BL = zeros(AbstractFloat, (Numθ₁, Numρ₂))
for nθ₁ = 1:length(θ₁)
	for nρ₂ = 1:length(ρ₂)
		BL[nθ₁, nρ₂] = OceanParameters.bottom_loss(θ₁[nθ₁]; ρ₂ = ρ₂[nρ₂])
	end
end

plot(rad2deg.(θ₁), BL)

cₛ = [6e2 4e2 2e2 0.0]
Numcₛ = length(cₛ)
BL = zeros(AbstractFloat, (Numθ₁, Numcₛ))
for nθ₁ = 1:length(θ₁)
	for ncₛ = 1:length(cₛ)
		BL[nθ₁, ncₛ] = OceanParameters.bottom_loss(θ₁[nθ₁]; cₛ = cₛ[ncₛ], αₚ = 0.)
	end
end

plot(rad2deg.(θ₁), BL)