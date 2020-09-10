module SonarEquations

using SpecialFunctions

"""
	SE = signal_excess_passive(SL, TL, NL, DI, DT)

Calculates the signal excess in a passive sonar scenario for source level `SL` (dB), transmission loss `TL` (dB), noise level `NL` (dB), directivity index `DI` (dB) and detection threshold `DT` (dB).
"""
function signal_excess_passive(SL, TL, NL, DI, DT)
	SE = SL - TL - NL + DI - DT
end

"""
	detection_index_gaussian(p_dtc, p_fal)

Calculates the detection index for a non-fluctuating signal in a Gaussian-distributed noise environment for a probability of detection `p_dtc` and probability of false alarm `p_fal`.
"""
function detection_index_gaussian(p_dtc, p_fal)
	d = 2(erfcinv(2p_dtc) - erfcinv(2p_fal))^2
end

function detection_index_gaussian(SL, TL, NL, B, t)
	d = B*t*((SL - TL)/(B*NL))^2
end

"""

"""
function detection_threshold(d, B, t)
	DT = 5log10(d*B/t)
end

"""
	probability_of_detection_gaussian(d, p_fal)

Calculates the probability of detection for a non-fluctuating signal in a Gaussian-distributed noise environment for detection index `d` and probability of false alarm `p_fal`.
"""
function probability_of_detection_gaussian(d, p_fal)
	p_dtc = erfc(erfcinv(2p_fal) - sqrt(d/2))/2
end

end
