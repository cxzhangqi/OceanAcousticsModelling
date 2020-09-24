## Preamble
using Plots

include("../AcousticPropagation.jl")

## n²-Linear Profile
c₀ = 1550
z₀ = 1e3
f = 2e3
src = AcousticPropagation.Entity(0, z₀, f)
ocn = AcousticPropagation.Medium((r, z) -> c₀/sqrt(1 + 2.4z/c₀), 3.5e3)
bty = AcousticPropagation.Boundary(z₀)
θ₀ = -acos(ocn.c(0, z₀)/ocn.c(0, 150))

# ray = AcousticPropagation.Ray(θ₀, src, ocn, bty)

# pt = plot(yaxis = :flip)
# plot!(ray.Sol, vars = (1, 2))
# display(pt)

δθ₀ = π/200
Nθ = 3
rays = AcousticPropagation.Rays(θ₀, δθ₀, Nθ, src, ocn, bty)

pt = plot(yaxis = :flip)
for nRay = 1:length(rays.rays)
	plot!(rays.rays[nRay].Sol, vars = (1, 2))
end
display(pt)