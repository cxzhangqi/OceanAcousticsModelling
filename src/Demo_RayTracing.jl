## Ray Tracing Demonstrations

## Preamble
using Plots
using StaticArrays

include("AcousticPropagation.jl")

## Demonstrative Scenario
# Altimetry
zAtiMin = -10
zAtiMax = 50
zAti(r) = zAtiMin + (zAtiMax - zAtiMin)*(sin(r/1e3) + 1.)/2

# Bathymetry
rPeak = 5e3
rMax = 10e3
zMax = 1e3
zMin = 8e2
Aᵣ = (2rPeak/3)^2/log((9zMax - 11zMin)/(10(zMax - zMin)))
zBty(r) = zMax - (zMax - zMin)*exp(-(r - rPeak)^2/4e5)

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

# Initial Conditions
r₀ = 0.
z₀ = (zBty(r₀) + zAti(r₀))/2
θ₀ = acos(c(r₀, z₀)/cMax).*(-1.5:0.5:1.5)

# Other
R = 2e4

RaySols = AcousticPropagation.helmholtz_eikonal.(θ₀, r₀, z₀, c, zAti, zBty, R)

r = range(0, R, length = 301)

pt = plot(
	legend = false,
	xaxis = "Range (m)",
	yaxis = ("Depth", :flip),
	title = "Rays")
plot!.(RaySols, vars = (1, 2),
	color = RGB(0, 0.5, 1))
plot!(r, zAti, color = :black)
plot!(r, zBty, color = :black)

display(pt)

savefig(pt, "RayTrace_FirstExample.png")

## Uniformly Increasing Celerity
zAtiVal = 0.
zBtyVal = 5e3
# zs = [0., 300., 1200., 2e3, 5e3]
# cs = [1520, 1500, 1515, 1495, 1545.]
# zs = [0, 5e3]
# cs = [1500, 1600]
# @show cMat = cat(zs, cs, dims = 2)
cFcn(r, z) = 1500 + 100z/5e3

r₀ = 0.0
z₀ = 20.0
θ₀ = deg2rad(10)
# θ₀ = range(atan(5000/25e3), atan(5000/50e3), length = 5)
# θ₀ = 5π/6
R = 250e3

RaySols = AcousticPropagation.helmholtz_eikonal.(θ₀, r₀, z₀, cFcn, zAtiVal, zBtyVal, R)

pt = plot(yaxis = :flip)
plot!.(RaySols, vars = (1, 2))
# scatter!(RaySols, vars = (1, 2))
display(pt)

##
zs = [0., 300., 1200., 2e3, 5e3]
cs = [1520, 1500, 1515, 1495, 1545.]

itp = interpolate((zs,), cs, Gridded(Linear()))
c(z) = itp(z)

z = range(0, 5e3, length = 101)
plot(z, c)

## Convergence Zone Propagation
zs = [0., 300., 1200., 2e3, 5e3]
cs = [1520, 1500, 1515, 1495, 1545.]
cMat = cat(zs, cs, dims = 2)

r₀ = 0.0
z₀ = 20.0
θ₀ = deg2rad(10)
R = 250e3

zAtiVal = 0
zBtyVal = 5e3

RaySol = AcousticPropagation.helmholtz_eikonal(θ₀, r₀, z₀, cMat, zAtiVal, zBtyVal, R)

pt = plot(yaxis = :flip)
plot!(RaySol, vars = (1, 2))
scatter!(RaySol, vars = (1, 2))
display(pt)

##
zs = [0., 300., 1200., 2e3, 5e3]
cs = [1520, 1500, 1515, 1495, 1545.]

itp = interpolate((zs,), cs, Gridded(Linear()))
c(z) = itp(z)

z = range(0, 5e3, length = 101)
plot(z, c)