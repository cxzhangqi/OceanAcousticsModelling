using Plots

include("SonarEquations.jl")

## ROC Curves: Gaussian-Gaussian
POD = range(0, 100, length = 100)
PFA = range(0, 100, length = 100)
contour(PFA, POD, (PFA, POD) -> SonarEquations.detection_index_gaussian(POD/100, PFA/100),
	fill = true,
	title = "Receiver Operating Characteristics Curves",
	xaxis = "POD (%)",
	yaxis = "PFA (%)")

## ROC Curves: Exponential-Exponential
