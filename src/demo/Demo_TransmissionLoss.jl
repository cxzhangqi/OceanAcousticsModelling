## Transmission Loss Demonstration

## Preamble
using Plots

## Lloyds Mirror
include("../LloydsMirror.jl")

SL = 50.
NL = 5.
B = 1.
t = 1.
DI = 20.
p_fal = 1e-2

c = 1500.0
f = 1e2
λ = c/f
r_src = 0.0
z_src = 2λ
r = range(0, 3e2, length = 1001)
z = range(0, 10λ, length = 501)
zTemp = 5λ

TL(r, z) = LloydsMirror.lloydsmirror_singlereflection.(c, f, r_src, r, z_src, z)

pt = heatmap(r, z, TL,
	seriescolor = cgrad(:jet, rev = true),
	title = "TL (dB)",
	yaxis = (:flip, "Depth (m)"),
	xaxis = ("Range (m)"))
display(pt)

savefig(pt, "img/SoundField_LloydsMirror_Simple.png")