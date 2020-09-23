## Parabolic Bathymetry
# Increase tolerance?
println("using Plots")
using Plots

println("Including AcousticPropagation.jl")
include("../AcousticPropagation.jl")

println("Defining Environment")
c = 250
R = 20e3
ocn = AcousticPropagation.Medium((r, z) -> c, R)
ati = AcousticPropagation.Boundary(r -> 0)
bty = AcousticPropagation.Boundary(r -> 2e-3*2.5e5sqrt(1 + r/c))
src = AcousticPropagation.Entity(0, 0)
θ₀ = π/4

println("Tracing Ray")
ray = AcousticPropagation.Ray(θ₀, src, ocn, bty, ati)

println("Plotting")
plot(yaxis = :flip)
plot!(range(0, R, length = 101), bty.z)
plot!(ray.Sol, vars = (1, 2))