"""
	AcousticPropagation

Calculates the sound pressure field for a cylindrical slice of the ocean.

TODO:
1. Provide multiple dispatch constructors for the range/depth dependent environment variables.
	a. Test.
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
* How to do tests for varying parameters?
"""
module AcousticPropagation

using Plots
using ForwardDiff
using DifferentialEquations
using LinearAlgebra

"""
	t_rfl::Vector = boundary_reflection(t_inc::Vector, t_bnd::Vector)

`t_inc`    tangent vector of incident ray
`t_bnd`    tangent vector of boundary
`t_rfl`    tangent vector of reflected ray
"""
function boundary_reflection(t_inc::Vector, t_bnd::Vector)
	n_bnd = [-t_bnd[2], t_bnd[1]]
# 	t_rfl = t_inc - 2(t_inc ⋅ n_bnd)*n_bnd
	t_rfl = t_inc - 2LinearAlgebra.dot(t_inc, n_bnd)*n_bnd
	return t_rfl
end

"""
	Boundary(z)

`z(r)`    depth of boundary in terms of range `r`
"""
struct Boundary
	z::Function
	dzdr::Function
	condition::Function
	affect!::Function
	function Boundary(z)
		dzdr(r) = ForwardDiff.derivative(z, r)
		condition(u, t, ray) = z(u[1]) - u[2]
		affect!(ray) = ray.u[3], ray.u[4] = boundary_reflection([ray.u[3], ray.u[4]], [1, dzdr(ray.u[1])])
		return new(z, dzdr, condition, affect!)
	end
end

"""
	Medium(c, R)

`c(r, z)`    sound speed at range `r` and depth `z`
`R`          maximum range of slice
"""
struct Medium
	c::Function
	∂c∂r::Function
	∂c∂z::Function
	R::Real
	function Medium(c::Function, R::Real)
		c_(x) = c(x[1], x[2])
		∇c_(x) = ForwardDiff.gradient(c_, x)
		∇c(r, z) = ∇c_([r, z])
		∂c∂r(r, z) = ∇c(r, z)[1]
		∂c∂z(r, z) = ∇c(r, z)[2]
		return new(c, ∂c∂r, ∂c∂z, R)
	end
end

"""
	Entity(r, z, θ)

`r`    range of entity
`z`    depth of entity
`θ`    angle of ray
"""
struct Entity
	r::Real
	z::Real
end

function acoustic_propagation_problem(
	θ₀::Real,
	Src::Entity,
	Ocn::Medium,
	Bty::Boundary,
	Ati::Boundary)

	function eikonal!(du, u, p, s)
		r = u[1]
		z = u[2]
		ξ = u[3]
		ζ = u[4]
		τ = u[5]
		du[1] = dr = Ocn.c(r, z)*ξ
		du[2] = dz = Ocn.c(r, z)*ζ
		du[3] = dξ = -Ocn.∂c∂r(r, z)/Ocn.c(r, z)^2
		du[4] = dζ = -Ocn.∂c∂z(r, z)/Ocn.c(r, z)^2
		du[5] = dτ = 1/Ocn.c(r, z)
	end

	rng_condition(u, t, ray) = Ocn.R/2 - abs(u[1] - Ocn.R/2)
	rng_affect!(ray) = terminate!(ray)
	CbRng = ContinuousCallback(rng_condition, rng_affect!)
	CbBty = ContinuousCallback(Bty.condition, Bty.affect!)
	CbAti = ContinuousCallback(Ati.condition, Ati.affect!)
	CbBnd = CallbackSet(CbRng, CbBty, CbAti)

	r₀ = Src.r
	z₀ = Src.z
	ξ₀ = cos(θ₀)/Ocn.c(r₀, z₀)
	ζ₀ = sin(θ₀)/Ocn.c(r₀, z₀)
	τ₀ = 0.
	u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀]

	TLmax = 100
	S = 10^(TLmax/10)
	sSpan = (0., S)

	prob_eikonal = ODEProblem(eikonal!, u₀, sSpan)

	return prob_eikonal, CbBnd
end

function solve_acoustic_propagation(prob_eikonal, CbBnd)
	@time RaySol = solve(prob_eikonal, callback = CbBnd)
	return RaySol
end

Base.broadcastable(m::Entity) = Ref(m)
Base.broadcastable(m::Medium) = Ref(m)
Base.broadcastable(m::Boundary) = Ref(m)

struct Ray
	θ₀
	Sol
	function Ray(θ₀, Src::Entity, Ocn::Medium, Bty::Boundary, Ati::Boundary)
		Prob, CbBnd = acoustic_propagation_problem(θ₀, Src, Ocn, Bty, Ati)
		Sol = solve_acoustic_propagation(Prob, CbBnd)
		return new(θ₀, Sol)
	end
end

# Old code below, to be updated.
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