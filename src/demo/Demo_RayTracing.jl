## Ray Tracing Demonstrations

## Preamble
using Plots

include("../AcousticPropagation.jl")

## Demonstrative Scenario
# Altimetry
zAtiMin = -10
zAtiMax = 50
zAti(r) = zAtiMin + (zAtiMax - zAtiMin)*(sin(r/1e3) + 1.)/2
ati = AcousticPropagation.Boundary(zAti)

# Bathymetry
rPeak = 5e3
rMax = 10e3
zMax = 1e3
zMin = 8e2
Aᵣ = (2rPeak/3)^2/log((9zMax - 11zMin)/(10(zMax - zMin)))
zBty(r) = zMax - (zMax - zMin)*exp(-(r - rPeak)^2/4e5)
bty = AcousticPropagation.Boundary(zBty)

# Ocean
cMin = 1500
cMax = 1600
C(r) = [1 zAti(r) zAti(r)^2
	1 (zAti(r) + zBty(r))/2 ((zAti(r) + zBty(r))/2)^2
	1 zBty(r) zBty(r)^2]
c_(r) = C(r)\[cMax, cMin, cMax]
c₀(r) = c_(r)[1]
c₁(r) = c_(r)[2]
c₂(r) = c_(r)[3]
c(r, z) = c₀(r) + c₁(r)*z + c₂(r)*z^2
ocn = AcousticPropagation.Medium(c, rMax)

# Initial Conditions
r₀ = 0.
z₀ = (zBty(r₀) + zAti(r₀))/2
src = AcousticPropagation.Entity.(r₀, z₀)

θ₀ = acos(c(r₀, z₀)/cMax).*(-1.5:0.125:1.5)
rays = AcousticPropagation.Ray.(θ₀, src, ocn, bty, ati)

println("Now plotting.")
pt = plot(yaxis = :flip)
r = range(0, rMax, length = 100)
# z = range(zAtiMin, zMax, length = 51)
plot!(r, zAti, label = "Altimetry")
plot!(r, zBty, label = "Bathymetry")
for nRay = 1:length(rays)
	plot!(rays[nRay].Sol, vars = (1, 2), label = "")
end
display(pt)

savefig(pt, "img/RayTrace_FirstExample.png")

## Uniformly Increasing Celerity
# Investigate difference between the interpolation and the function
ati = AcousticPropagation.Boundary(0)
bty = AcousticPropagation.Boundary(5e3)

# zs = [0., 300., 1200., 2e3, 5e3]
# cs = [1520, 1500, 1515, 1495, 1545.]

zs = [0, 5e3]
cs = [1500, 1600]

R = 250e3
cMat = vcat([0 0 R], hcat(zs, cs, cs))
ocnVec = AcousticPropagation.Medium(cMat, R)

c(r, z) = 1500 + 100(z - zAtiVal)/(zBtyVal - zAtiVal)
ocnFcn = AcousticPropagation.Medium(c, R)

h = heatmap(range(0, R, length = 51),
	range(zs[1], zs[end], length = 101),
	(r, z) -> ocnFcn.c(r, z) - ocnVec.c(r, z))
display(h)

src = AcousticPropagation.Entity(0., 20.)
θ₀ = deg2rad(10)
# θ₀ = range(atan(5000/25e3), atan(5000/50e3), length = 5)
# θ₀ = 5π/6

ray = AcousticPropagation.Ray(θ₀, src, ocn, bty, ati)

pt = plot(yaxis = :flip)
plot!(ray.Sol, vars = (1, 2))
display(pt)

## Parabolic Bathymetry
# SO:
# * Parabolic works for angle-controlled reflection.
# * General works for vector-controlled reflection.
cVal = 250
R = 20e3
ocn = AcousticPropagation.Medium((r, z) -> cVal, R)
ati = AcousticPropagation.Boundary(r -> 0)
bty = AcousticPropagation.Boundary(r -> 2e-3*2.5e5sqrt(1 + r/cVal))
src = AcousticPropagation.Entity(0, 0)
θ₀ = range(atan(5e3/2e3), atan(5e3/20e3), length = 10)

rays = AcousticPropagation.Ray.(θ₀, src, ocn, bty, ati)

pt = plot(yaxis = :flip)
plot!(range(0, R, length = 101), bty.z)
# plot!(ray.Sol, vars = (1, 2))
for nRay = 1:length(rays)
	plot!(rays[nRay].Sol, vars = (1, 2), label = "")
end
display(pt)