##
using LinearAlgebra
include("AcousticPropagation.jl")
# Altimetry
zAtiMin = -10
zAtiMax = 50
zAti(r) = zAtiMin + (zAtiMax - zAtiMin)*(sin(r/1e3) + 1.)/2

# Bathymetry
rPeak = 5e3
rMax = 10e3
zMax = 1e3
zMin = 7e2
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
S = 2e4 # Figure out how to replace with condition

RaySols = helmholtz_eikonal_transport.(θ₀, r₀, z₀, c, zAti, zBty, S)

using Plots

pt = plot(yaxis = :flip)
plot!.(RaySols, vars = (1, 2))
display(pt)


##
using Images
using FileIO
SvpImg = load("img/SVP_ConvergenceZonePropagation.png")

using Images
using TestImages
using FileIO
img = testimage("mandrill")

##
using Plots; gr()
using SpecialFunctions
using ColorSchemes

B = 1.
t = 1.
SL = 52.
NL = 5.
p_fal = 1e-4

c = 1500.
f = 100.
λ = c/f
k = 2π/λ

r₀ = 0.
z₀ = 100.

R(r, z) = sqrt.((r - r₀)^2 + (z - z₀)^2)
p(r, z) = R(r, z)*exp(im*k*R(r, z))
TL(r, z) = 10log10(abs(p(r, z)))

d(r, z) = B*t*((SL - TL(r, z))/(B*NL))^2
p_dtc(r, z) = erfc(erfcinv(2p_fal) - sqrt(d(r, z)/2))/2
POD(r, z) = 100p_dtc(r, z)

r = range(1., 1e3, length = 100)
z = range(1., 400., length = 50)
contour(r, z, POD,
	fill = true,
	seriescolor = :jet,
	xaxis = "Range (m)",
	yaxis = ("Depth (m)", :flip))