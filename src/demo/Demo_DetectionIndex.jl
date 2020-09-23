using Plots

include("../SonarEquations.jl")

## ROC Curves: Gaussian-Gaussian
POD = range(0, 100, length = 101)
PFA = 10.0.^range(-4, 1, length = 101)
d(p_dtc, p_fal) = SonarEquations.detection_index_gaussian(p_dtc, p_fal)
d_PC_dB(PFA, POD) = 10log10(d(POD/100, PFA/100))
pt = contour(PFA, POD, d_PC_dB,
	title = "ROC Curves: Gaussian",
	xaxis = ("PFA (%)", :log10),
	yaxis = "POD (%)",
	colorbar_title = "Detection Index",
	levels = [0, 3, 6, 9, 10, 12, 15])

display(pt)

savefig(pt, "img/DetectionIndex_Gaussian.png")

# plot(POD/100, p_dtc -> d(p_dtc, 1e-1))
# plot(PFA/100, p_fal -> d(0.5, p_fal), xaxis = :log10)
	
## ROC Curves: Exponential-Exponential
