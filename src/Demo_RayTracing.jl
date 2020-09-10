using Plots

include("AcousticPropagation.jl")

## Convergenze-Zone Propagation (soon)
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

RaySols = AcousticPropagation.helmholtz_eikonal_transport.(θ₀, r₀, z₀, c, zAti, zBty, S)

pt = plot(yaxis = :flip)
plot!.(RaySols, vars = (1, 2))
display(pt)

## Convergenze-Zone Propagation
function LineFcn(x₁, x₂, y₁, y₂)
	y(x) = y₁ + (y₂ - y₁)/(x₂ - x₁)*(x - x₁)
end

zAti(r) = 0.
zBty(r) = 5e3
function c(r, z)
	zs = [0., 300., 1200., 2e3, 5000.]
	cs = [1520, 1500, 1515, 1495, 1545.]
	
	n₋ = findlast(zs .≤ z)
	n₊ = findfirst(z .< zs)

	cMin = cs[n₋]
	cMax = cs[n₊]
	zMin = zs[n₋]
	zMax = zs[n₊]
	cFcn(z) = cMin + (cMax - cMin)*(z - zMin)*(zMax - zMin)

	return cFcn(z)
end
r₀ = 0.
z₀ = 20.
θ₀ = atan(1/6)
S = 140e3

RaySols = AcousticPropagation.helmholtz_eikonal_transport(θ₀, r₀, z₀, c, zAti, zBty, S)

pt = plot(yaxis = :flip)
plot!(RaySols, vars = (1, 2))
display(pt)
