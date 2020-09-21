using Plots

include("SonarEquations.jl")

## ROC Curves: Gaussian-Gaussian
POD = range(0, 100, length = 100)
PFA = range(0, 100, length = 100)
pt = contour(PFA, POD, (PFA, POD) -> SonarEquations.detection_index_gaussian(POD/100, PFA/100),
	title = "Receiver Operating Characteristics Curves",
	xaxis = "POD (%)",
	yaxis = "PFA (%)",
	colorbar_title = "Detection Index")

savefig(pt, "img/DetectionIndex_Gaussian.png")
	
## ROC Curves: Exponential-Exponential
