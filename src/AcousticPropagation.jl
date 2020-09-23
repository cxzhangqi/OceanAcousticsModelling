"""
	AcousticPropagation

Calculates the sound pressure field for a cylindrical slice of the ocean.

TODO:
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
using ForwardDiff
using DifferentialEquations
using LinearAlgebra
using Interpolations

"""

"""
function InterpolatingFunction(rng, val)
	Itp = LinearInterpolation(rng, val, extrapolation_bc = Flat())
	return ItpFcn(r) = Itp(r)
end

"""

"""
function InterpolatingFunction(rng, dpt, val)
	display(rng)
	display(dpt)
	Itp = LinearInterpolation((dpt, rng), val, extrapolation_bc = Flat())
	return ItpFcn(r, z) = Itp(z, r)
end

"""
	t_rfl::Vector = boundary_reflection(t_inc::Vector, t_bnd::Vector)

`t_inc`    tangent vector of incident ray
`t_bnd`    tangent vector of boundary
`t_rfl`    tangent vector of reflected ray

TODO:
* Reconcile the two versions via error checking.
"""
function boundary_reflection(t_inc::Vector, t_bnd::Vector)
	# works for parabolic boundary
	MyAngle(tng) = atan(tng[2]/tng[1])
	θ_inc = MyAngle(t_inc)
	θ_bnd = MyAngle(t_bnd)

	c = cos(θ_inc)/t_inc[1]

	θ_inc_flat = θ_inc - θ_bnd
	θ_rfl_flat = -θ_inc_flat
	θ_rfl = θ_rfl_flat + θ_bnd

	return [cos(θ_rfl), sin(θ_rfl)]/c
end
# function boundary_reflection(t_inc::Vector, t_bnd::Vector)
# 	# works for generic boundary
# 	n_bnd = [-t_bnd[2], t_bnd[1]]
# # 	t_rfl = t_inc - 2(t_inc ⋅ n_bnd)*n_bnd
# 	t_rfl = t_inc - 2LinearAlgebra.dot(t_inc, n_bnd)*n_bnd

# 	MyAngle(tng) = atand(tng[2]/tng[1])
# 	θ_inc = MyAngle(t_inc)
# 	θ_bnd = MyAngle(t_bnd)
# 	θ_rfl = MyAngle(t_rfl)

# 	# if !isapprox(θ_inc - θ_bnd, θ_bnd - θ_rfl; atol = 1e-6)
# 	# 	println(t_inc)
# 	# 	println(θ_inc)
# 	# 	println(t_bnd)
# 	# 	println(θ_bnd)
# 	# 	println(t_rfl)
# 	# 	println(θ_rfl)
# 	# 	println(θ_inc - θ_bnd)
# 	# 	println(θ_bnd - θ_rfl)
# 	# 	error("Angle calculated incorrectly.")
# 	# end

# 	return t_rfl
# end

"""
	Boundary(z)

`z(r)`    depth of boundary in terms of range `r`
"""
struct Boundary
	z::Function
	dzdr::Function
	condition::Function
	affect!::Function
end
function Boundary(z::Function)
	dzdr(r) = ForwardDiff.derivative(z, r)
	condition(u, t, ray) = z(u[1]) - u[2]
	affect!(ray) = ray.u[3], ray.u[4] = boundary_reflection([ray.u[3], ray.u[4]], [1, dzdr(ray.u[1])])
	return Boundary(z, dzdr, condition, affect!)
end
function Boundary(rz::AbstractArray)
	r_ = [rng for rng ∈ rz[:, 1]]
	z_ = [dpt for dpt ∈ rz[:, 2]]
	zFcn = InterpolatingFunction(r_, z_)
	return Boundary(zFcn)
end
function Boundary(z::Real)
	zFcn(r) = z
	return Boundary(zFcn)
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
	∂²c∂r²::Function
	∂²c∂r∂z::Function
	∂²c∂z²::Function
	R::Real
end
function Medium(c::Function, R::Real)
	c_(x) = c(x[1], x[2])
	∇c_(x) = ForwardDiff.gradient(c_, x)
	∇c(r, z) = ∇c_([r, z])
	∂c∂r(r, z) = ∇c(r, z)[1]
	∂c∂z(r, z) = ∇c(r, z)[2]

	∂c∂r_(x) = ∂c∂r(x[1], x[2])
	∇∂c∂r_(x) = ForwardDiff.gradient(∂c∂r_, x)
	∇∂c∂r(r, z) = ∇∂c∂r_([r, z])

	∂c∂z_(x) = ∂c∂z(x[1], x[2])
	∇∂c∂z_(x) = ForwardDiff.gradient(∂c∂r_, x)
	∇∂c∂z(r, z) = ∇∂c∂z_([r, z])

	∂²c∂r²(r, z) = ∇∂c∂r(r, z)[1]
	∂²c∂r∂z(r, z) = ∇∂c∂r(r, z)[2]
	∂²c∂z²(r, z) = ∇∂c∂z(r, z)[2]

	return Medium(c, ∂c∂r, ∂c∂z, ∂²c∂r², ∂²c∂r∂z, ∂²c∂z², R)
end
function Medium(c::AbstractArray, R::Real)
	r_ = [rc for rc ∈ c[1, 2:end]]
	z_ = [zc for zc ∈ c[2:end, 1]]
	c_ = c[2:end, 2:end]
	
	cFcn = InterpolatingFunction(r_, z_, c_)
	return Medium(cFcn, R)
end
function Medium(c::Real, R::Real)
	cFcn(r, z) = c
	return Medium(cFcn, R)
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
	q = u[6]
	p = u[7]

	∂²c∂n²(r, z) = Ocn.c(r, z)^2*(
		Ocn.∂²c∂r²(r, z)*ζ^2
		- 2Ocn.∂²c∂r∂z(r, z)*ξ*ζ
		+ Ocn.∂²c∂z²(r, z)*ξ^2
	)

	du[1] = drds = Ocn.c(r, z)*ξ
	du[2] = dzds = Ocn.c(r, z)*ζ
	du[3] = dξds = -Ocn.∂c∂r(r, z)/Ocn.c(r, z)^2
	du[4] = dζds = -Ocn.∂c∂z(r, z)/Ocn.c(r, z)^2
	du[5] = dτds = 1/Ocn.c(r, z)
	du[6] = dqds = Ocn.c(r, z)*p
	du[7] = dpds = ∂²c∂n²(r, z)/Ocn.c(r, z)^2*q
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
τ₀ = 0.0
q₀ = 0.0
p₀ = 1.0/Ocn.c(r₀, z₀)
u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, q₀, p₀]

TLmax = 100
S = 10^(TLmax/10)
sSpan = (0., S)

prob_eikonal = ODEProblem(eikonal!, u₀, sSpan)

return prob_eikonal, CbBnd
end

function solve_acoustic_propagation(prob_eikonal, CbBnd)
	@time RaySol = solve(prob_eikonal, callback = CbBnd, reltol=1e-8, abstol=1e-8)
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

end