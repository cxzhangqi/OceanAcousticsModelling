## Demo: Sonar Equations
# Demonstrate calculation of the sonar equations, with focus on signal excess and probability of detection.
using Plots; gr()
using SpecialFunctions
using ColorSchemes

include("SonarEquations.jl")

## Simple Scenario
B = 1. # bandwidth
t = 1. # integration time
SL = 52. # source level
NL = 5. # noise level
p_fal = 1e-4 # probability of false alarm

c = 1500. # sound speed
f = 100. # frequency
λ = c/f # wavelength
k = 2π/λ # wavenumber

r₀ = 0. # source range
z₀ = 100. # source depth

R(r, z) = sqrt.((r - r₀)^2 + (z - z₀)^2) # distance from source
p(r, z) = R(r, z)*exp(im*k*R(r, z)) # complex pressure
TL(r, z) = 10log10(abs(p(r, z))) # transmission loss
DI = 20 # directivity index
# DT(r, z) = 5log10(d(r, z)*B/t)
DT(r, z) = SonarEquations.detection_threshold(d(r, z), B, t)
# SE(r, z) = SL - TL(r, z) - NL + DI - DT(r, z)
SE(r, z) = SonarEquations.signal_excess_passive(SL, TL(r, z), NL, DI, DT(r, z))

# d(r, z) = B*t*((SL - TL(r, z))/(B*NL))^2
d(r, z) = SonarEquations.detection_index_gaussian(SL, TL(r, z), NL, B, t)
# p_dtc(r, z) = erfc(erfcinv(2p_fal) - sqrt(d(r, z)/2))/2
# POD(r, z) = 100p_dtc(r, z)
POD(r, z) = 100SonarEquations.probability_of_detection_gaussian(d(r, z), p_fal)

r = range(1., 1e3, length = 100)
z = range(1., 400., length = 50)

pd = contour(r, z, d,
	fill = true,
	seriescolor = :jet,
	title = "d",
	yaxis = ("Depth (m)", :flip))

pTL = contour(r, z, TL,
	fill = true,
	seriescolor = cgrad(:jet, rev = true),
	title = "TL (dB)",
	yaxis = :flip)

pSE = contour(r, z, SE,
	fill = true,
	seriescolor = :jet,
	title = "SE (dB)",
	xaxis = "Range (m)",
	yaxis = ("Depth", :flip))

pPOD = contour(r, z, POD,
	fill = true,
	seriescolor = :jet,
	title = "POD (%)",
	xaxis = "Range (m)",
	yaxis = :flip)

l = @layout [a b; c d]

pt = plot(pd, pTL, pSE, pPOD,
	layout = l)

savefig(pt, "img/SonarEqs_SimplePropagation.png")

## Simple Lloyd's Mirror

## Multipath Propagation

