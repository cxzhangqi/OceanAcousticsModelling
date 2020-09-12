"""
	AcousticPropagation

Calculates the sound pressure field for a cylindrical slice of the ocean.

TODO:
1. Divide into `problem` definition and `solving` definition for interaction with intermediate data, e.g. own plots of slice.
	0. Copy.
	a. Produce struct with multiple dispatch constructors to define environmental features.
	b. Divide solving function.
	c. Test on working
2. Define mature plotting function.
	a. Change output type to own struct `Ray`.
	b. Extend base plotting function to plot with `Ray` input?
3. Test ray tracing on scenarios.
	a. Change computation condition to range instead of arc length.
	b. Create mature plotting function.
	c. Convert to structs? Maybe later, for 2.0.
4. Implement amplitude calculation.
5. Produce pressure field output.
6. Combine with sonar equations.
7. Test pressure and TL field on scenarios.
8. Create as formal Julia package 1.0 (or 0.1?) as proof of concept.

Questions:
* How to rename library
* LinearAlgebra: dot doesn't work: `UndefVarError: ⋅ not defined`
* How to do tests for varying parameters?
"""
module AcousticPropagation

# Dependencies
using Interpolations:LinearInterpolation
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

Calculates the ray path geometry, stored in `RaySol` for initial grazing angle `θ₀` (rad), initial range `r₀`, initial depth `z₀`, sound ocean speed field `c`, ocean surface altimetry `zAti`, ocean bottom bathymetry `zBty` and maximum range `R`.

TODO:
* Consider zero-range reflection for cylindrically symmetric environment
* Calculate amplitude
* Calculate pressure
"""
function helmholtz_eikonal(θ₀::Real, r₀::Real, z₀::Real, c::Function, zAti::Function, zBty::Function, R::Real)
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

	τ₀ = 0.
	ξ₀ = cos(θ₀)/c(r₀, z₀)
	ζ₀ = sin(θ₀)/c(r₀, z₀)
	u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀]
	TLmax = 100
	S = 10^(TLmax/10)
	tspan = (0., S)
	prob = ODEProblem(eikonal!, u₀, tspan)
	@time RaySol = solve(prob, 
		callback = cb_bnd)
	return RaySol
end
function helmholtz_eikonal(θ₀::Real, r₀::Real, z₀::Real, c, zAti, zBty, R::Real)
	cFcn = parse_speed(c, R)
	zAtiFcn = parse_boundary(zAti, R)
	zBtyFcn = parse_boundary(zBty, R)
	return helmholtz_eikonal(θ₀, r₀, z₀, cFcn, zAtiFcn, zBtyFcn, R)
end
# function helmholtz_eikonal(θ₀::Real, r₀::Real, z₀::Real, c::AbstractArray, zAti::AbstractArray, zBty::AbstractArray, R::Real)

# 	itpBty = interpolate((zBty[:, 1],), zBty[:, 2], Gridded(Linear()))
# 	zBtyFcn(r) = itpBty(r)

# 	itpAti = interpolate((zAti[:, 1],), zAti[:, 2], Gridded(Linear()))
# 	zAtiFcn(r) = itpAti(r)

# 	return helmholtz_eikonal(θ₀, r₀, z₀, cFcn, zAtiFcn, zBtyFcn, R)
# end
# function helmholtz_eikonal(θ₀::Real, r₀::Real, z₀::Real, c::Real, zAti::Real, zBty::Real, R::Real)
# 	println("Used dispatch!")
# 	cFcn(r, z) = c
# 	zAtiFcn(r) = zAti
# 	zBtyFcn(r) = zBty
# 	return helmholtz_eikonal(θ₀, r₀, z₀, cFcn, zAtiFcn, zBtyFcn, R)
# end

"""
	parse_speed(c::AbstractArray)

Available input types:
* TODO
"""
function parse_speed(c::Function, R::Real)
	return c
end
function parse_speed(c::AbstractArray, R::Real)
	# Maybe use first element to indicate interpolation degree?
	if size(c, 2) == 2
		println("=== Interpolate 2 ===")
		dep = c[:, 1]
		pushfirst!(dep, -R)
		push!(dep, 1.5R)
		spd = c[:, 2]
		pushfirst!(spd, spd[1])
		push!(spd, spd[end])
		itpSpd = LinearInterpolation((dep,), spd)
		return cFcn(r, z) = itpSpd(z)
		# cFcnTemp(z) = itpSpd(z)
		# range = 
		# plot()
	elseif size(c, 2) > 2
		# TODO:
		# * Implement as vector of arrays.
		# * First element is 
	end
end
function parse_speed(c::Real, R::Real)
	return cFcn(r, z) = c
end

"""
	parse_boundary(zBnd)

Available input types:
* TODO
"""
function parse_boundary(zBnd::Function, R::Real)
	return zBnd
end
function parse_boundary(zBnd::AbstractArray, R::Real)
	rng = zBnd[:, 1]
	pushfirst!(rng, -R)
	push!(rng, 1.5R)
	bnd = zBnd[:, 2]
	pushfirst!(bnd, bnd[1])
	push!(bnd, bnd[end])
	itpBnd = interpolate((rng,), bnd, Gridded(Linear()))
	return zBndFcn(r) = itpBnd(r)
end
function parse_boundary(zBnd::Real, R::Real)
	return zBndFcn(r) = zBnd
end

end