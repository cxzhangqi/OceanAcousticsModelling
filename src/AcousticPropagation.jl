"""
	AcousticPropagation

Calculates the sound field for a cylindrical slice of the ocean.

Questions:
* How to rename library
* LinearAlgebra: dot doesn't work: `UndefVarError: ⋅ not defined`
* How to do tests for varying parameters?
"""
module AcousticPropagation

# Dependencies
using LinearAlgebra: dot
using ForwardDiff: ForwardDiff
using DifferentialEquations:
ContinuousCallback,
CallbackSet,
ODEProblem,
solve
# using LinearAlgebra
# using ForwardDiff
# using DifferentialEquations

"""
	t_rfl = boundary_reflection(t_inc::Vector, t_bnd::Vector)

Calculates the reflection vector `t_rfl` for a incident vector `t_inc` reflecting off the boundary vector `t_bnd`.
"""
function boundary_reflection(t_inc::Vector, t_bnd::Vector)
	n_bnd = [-t_bnd[2], t_bnd[1]]
	@show t_inc
	@show typeof(t_inc)
	@show n_bnd
	@show typeof(n_bnd)
	t_rfl = t_inc - 2(t_inc ⋅ n_bnd)*n_bnd
	return t_rfl
end

"""
	RaySol = helmholtz_eikonal_transport(θ₀, r₀, z₀, c, zAti, zBty, S)

Calculates the ray path geometry, stored in `RaySol` for initial grazing angle `θ₀`, initial range `r₀`, initial depth `z₀`, sound ocean speed field `c`, ocean surface altimetry `zAti`, ocean bottom bathymetry `zBty` for a path of length `S`.

TODO:
* Calculate amplitude
* Calculate pressure
"""
function helmholtz_eikonal_transport(θ₀, r₀, z₀, c, zAti, zBty, S)
	println(θ₀)
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
		dr = c(r, z)*ξ
		dz = c(r, z)*ζ
		dξ = -∂c∂r(r, z)/c(r, z)^2
		dζ = -∂c∂z(r, z)/c(r, z)^2
		dτ = 1/c(r, z)
		du[1] = dr
		du[2] = dz
		du[3] = dξ
		du[4] = dζ
		du[5] = dτ
	end

	dzdrAti(r) = ForwardDiff.derivative(zAti, r)
	ati_condition(u, t, ray) = zAti(u[1]) - u[2]
	function ati_affect!(ray)
		r = ray.u[1]
		z = ray.u[2]
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
		z = ray.u[2]
		ξᵢ = ray.u[3]
		ζᵢ = ray.u[4]
		ξᵣ, ζᵣ = boundary_reflection([ξᵢ, ζᵢ], [1, dzdrBty(r)])
		ray.u[3] = ξᵣ
		ray.u[4] = ζᵣ
	end
	cb_bty = ContinuousCallback(bty_condition, bty_affect!)

	cb_bnd = CallbackSet(cb_ati, cb_bty)

	τ₀ = 0
	ξ₀ = cos(θ₀)/c(r₀, z₀)
	ζ₀ = sin(θ₀)/c(r₀, z₀)
	u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀]
	tspan = (0., S)
	prob = ODEProblem(eikonal!, u₀, tspan)
	@time RaySol = solve(prob, callback = cb_bnd)
	return RaySol
end

end