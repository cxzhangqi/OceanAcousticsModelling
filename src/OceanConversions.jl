"""
	x	range
	z	depth (+ve downards)
	p	pressure
"""
module OceanConversions
export PressureToTransmissionLoss

function PressureToTransmissionLoss(p::Complex)
	TL = PressureToTransmissionLoss(abs(p))
end

function PressureToTransmissionLoss(p::Real)
	I = p^2
	TL = IntensityToTransmissionLoss(I)
end

function IntensityToTransmissionLoss(I::Real)
	TL = -10log10(I)
end

function TransmissionLoss_SphericalSpreading(r)
	TL = 20log10(r)
end