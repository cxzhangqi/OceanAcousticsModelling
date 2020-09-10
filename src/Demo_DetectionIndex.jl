using Plots

include("SonarEquations.jl")

p_dtc = 10 .^range(-2, log10(0.99), length = 100)
p_fal = 10 .^range(-6, log10(0.99), length = 100)
# p_dtc = log10.(range(10^(1e-2), 10^(0.99), length = 100))
# p_fal = log10.(range(10^(1e-6), 10^(1e-1), length = 100))
contour(p_dtc, p_fal, SonarEquations.detection_index_gaussian,
	fill = true)

