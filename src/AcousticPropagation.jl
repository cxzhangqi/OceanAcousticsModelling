"""
	AcousticPropagation

Calculates the sound pressure field for a cylindrical slice of the ocean.

TODO:
1. Test ray tracing on scenarios.
	a. Change computation condition to range instead of arc length.
	b. Create mature plotting function.
	c. Convert to structs? Maybe later, for 2.0.
2. Implement amplitude calculation.
3. Produce pressure field output.
4. Combine with sonar equations.
5. Test pressure and TL field on scenarios.
6. Create as formal Julia package 1.0 (or 0.1?) as proof of concept.

Questions:
* How to rename library
* LinearAlgebra: dot doesn't work: `UndefVarError: ⋅ not defined`
* How to do tests for varying parameters?
"""
module AcousticPropagation

# Dependencies
using StaticArrays
# using LinearAlgebra: dot
using LinearAlgebra
using ForwardDiff: ForwardDiff
using DifferentialEquations:
ContinuousCallback,
CallbackSet,
ODEProblem,
solve,
terminate!
# using LinearAlgebra
# using ForwardDiff
# using DifferentialEquations

"""
	t_rfl = boundary_reflection(t_inc::Vector, t_bnd::Vector)

Calculates the reflection vector `t_rfl` for a incident vector `t_inc` reflecting off the boundary vector `t_bnd`.
"""
function boundary_reflection(t_inc::Vector, t_bnd::Vector)
	n_bnd = [-t_bnd[2], t_bnd[1]]
	t_rfl = t_inc - 2(t_inc ⋅ n_bnd)*n_bnd
	return t_rfl
end

"""
	RaySol = helmholtz_eikonal_transport(θ₀, r₀, z₀, c, zAti, zBty, S)

Calculates the ray path geometry, stored in `RaySol` for initial grazing angle `θ₀`, initial range `r₀`, initial depth `z₀`, sound ocean speed field `c`, ocean surface altimetry `zAti`, ocean bottom bathymetry `zBty` for a path of length `S` (to be deprecated) and maximum range `R`.

TODO:
* Replace stopping condition
* Calculate amplitude
* Calculate pressure
"""
function helmholtz_eikonal(θ₀, r₀, z₀, c, zAti, zBty, S, R)
	c_(x) = c(x[1], x[2])
	∇c_(x) = ForwardDiff.gradient(c_, x)
	∇c(r, z) = ∇c_([r, z])
	∂c∂r(r, z) = ∇c(r, z)[1]
	∂c∂z(r, z) = ∇c(r, z)[2]

	"""
		eikonal!(du, u, p, s)
	
	Calculates the change of parameters `du` of `u` for the eikonal equation.

	TODO:
	* Check if `p` and `s` need to be input variables in the function definition.
	"""
	function eikonal!(du, u, p, s)
		r = u[1]
		z = u[2]
		ξ = u[3]
		ζ = u[4]
		τ = u[5]
		du[1] = dr = c(r, z)*ξ
		du[2] = dz = c(r, z)*ζ
		du[3] = dξ = -∂c∂r(r, z)/c(r, z)^2
		du[4] = dζ = -∂c∂z(r, z)/c(r, z)^2
		du[5] = dτ = 1/c(r, z)
	end

	dzdrAti(r) = ForwardDiff.derivative(zAti, r)
	ati_condition(u, t, ray) = zAti(u[1]) - u[2]
	function ati_affect!(ray)
		r = ray.u[1]
		# z = ray.u[2]
		ξᵢ = ray.u[3]
		ζᵢ = ray.u[4]
		ξᵣ, ζᵣ = boundary_reflection([ξᵢ, ζᵢ], [1, dzdrAti(r)])
		ray.u[3] = ξᵣ
		ray.u[4] = ζᵣ
	end
	cb_ati = ContinuousCallback(ati_condition, ati_affect!)

	dzdrBty(r) = ForwardDiff.derivative(zBty, r)
	bty_condition(u, t, ray) = u[2] - zBty(u[1])
	function bty_affect!(ray)
		r = ray.u[1]
		# z = ray.u[2]
		ξᵢ = ray.u[3]
		ζᵢ = ray.u[4]
		ξᵣ, ζᵣ = boundary_reflection([ξᵢ, ζᵢ], [1, dzdrBty(r)])
		ray.u[3] = ξᵣ
		ray.u[4] = ζᵣ
	end
	cb_bty = ContinuousCallback(bty_condition, bty_affect!)

	r_condition(u, t, ray) = R/2 - abs(u[1] - R/2)
	r_affect!(ray) = terminate!(ray)
	cb_range = ContinuousCallback(r_condition, r_affect!)

	cb_bnd = CallbackSet(cb_ati, cb_bty, cb_range)

	τ₀ = 0
	ξ₀ = cos(θ₀)/c(r₀, z₀)
	ζ₀ = sin(θ₀)/c(r₀, z₀)
	u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀]
	# TLmax = 100
	# S = 10^(TLmax/10)
	tspan = (0., S)
	prob = ODEProblem(eikonal!, u₀, tspan)
	@time RaySol = solve(prob, 
		callback = cb_bnd)
	return RaySol
end

end