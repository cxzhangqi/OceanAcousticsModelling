##
using DifferentialEquations, Plots, ForwardDiff

struct Boundary
	z
	dzr
	function Boundary(z)
		dzr(r) = ForwardDiff.derivative(z, r)
		return new(z, dzr)
	end
end

struct Entity
	r
	z
end

struct Medium
	c
	∂c_∂r
	∂c_∂z
	function Medium(c)
		c_(x) = c(x[1], x[2])
		∇c(x) = ForwardDiff.gradient(c_, x)
		∇c(r, z) = ∇c([r, z])
		∂c_∂r(r, z) = ∇c(r, z)[1]
		∂c_∂z(r, z) = ∇c(r, z)[2]
		return new(c, ∂c_∂r, ∂c_∂z)
	end
end

"""
Features:
* Range independence.
* Depth dependence for sound speed and bathymetry.

TODO:
* Add ray pressure calculation including solving the transport equation.
* Find out how to do automatic partial derivatives.
* Find out how to make solver stop upon reaching a range instead of time dependence.
* Check initial depth is between bathymetry and altimetry.
* Multiple dispatch on struct contructor methods.
* Report issue to DiffEqBase.
* Generalise inputs.
* Replicate Jensen ray tracing.
"""
function RayTrace(θ₀, Ocn, Src, Bty, Ati, T)
		ξ₀ = cos.(θ₀)/Ocn.c(Src.r, Src.z)
		ζ₀ = sin.(θ₀)/Ocn.c(Src.r, Src.z)
		u0 = [Src.r, Src.z, ξ₀, ζ₀]

		function BoundaryReflection!(Ray, Bnd::Boundary)
			@show θ_bnd = atan(Bnd.dzr(Ray.u[2]))
			@show θ_inc = asin(Ray.u[4]*Ocn.c(Src.r, Src.z))
			@show θ_out = 2θ_bnd - θ_inc
			Ray.u[4] = sin(θ_out)/Ocn.c(Src.r, Src.z)
			Ray.u[3] = cos(θ_out)/Ocn.c(Src.r, Src.z)
			@show rad2deg(θ_inc - θ_bnd)
			@show rad2deg(θ_out - θ_bnd)
		end

		# Altimetry
		Aty_condition(u, t, integrator) = u[2] - Ati.z
		Aty_affect!(integrator) = integrator.u[4] *= -1
		Aty_cb = ContinuousCallback(Aty_condition, Aty_affect!)
		
		# Bathymetry
		Bty_condition(u, t, integrator) = u[2] - Bty.z(u[1])
		# Bty_affect!(integrator) = integrator.u[4] *= -1
		Bty_affect!(integrator) = BoundaryReflection!(integrator, Bty)
		Bty_cb = ContinuousCallback(Bty_condition, Bty_affect!)
		
		Bnd_cb = CallbackSet(Aty_cb, Bty_cb)

		function eikonal!(du, u, p, t)
			r = u[1]
			z = u[2]
			ξ = u[3]
			ζ = u[4]
			du[1] = Ocn.c(r, z)*ξ
			du[2] = Ocn.c(r, z)*ζ
			du[3] = -Ocn.∂c_∂r(r, z)/Ocn.c(r, z)^2
			du[4] = -Ocn.∂c_∂z(r, z)/Ocn.c(r, z)^2
		end

		tspan = (0., T)
		prob = ODEProblem(eikonal!, u0, tspan)

		sol = solve(prob, callback = Bnd_cb)
end

T = 1e4
cMin = 1.5e3
cMax = 1.6e3
zMax = 1e3
cz = [0 0 1; 2.5e5 5e2 1; 1e6 1e3 1.0]\[cMax; cMin; cMax]
# cOcn(r, z) = cMin + ((cz[1]*z^2 + cz[2]*z + cz[3]) - cMin)
cOcn(r, z) = cMin + (cMax - cMin)*abs(z - zMax/2)/(zMax/2)
r₀ = 0.
z₀ = 250.
zBeg = 1e3
zEnd = 7.5e2
# zBty(r) = zBeg + (zEnd - zBeg)*r/10e3
zBty(r) = 1e3 - 2e2exp(-(r - 3e3)^2/1e6)

Bty = Boundary(zBty)
Ati = Boundary(0.0)
Src = Entity(r₀, z₀)
Ocn = Medium(cOcn)
θ₀_max = acos(cOcn(r₀, z₀)/cMax)
# θ₀_vec = (-1.5:0.25:1.5).*θ₀_max
θ₀_vec = θ₀_max

let pt, rng
for n = 1:length(θ₀_vec)
	RayPath = RayTrace(θ₀_vec[n], Ocn, Src, Bty, Ati, T)
	if n == 1
		rng = RayPath[1,:]
		pt = plot(RayPath,
			vars = (1,2), 
			xaxis = "Range (m)",
			yaxis = ("Depth (m)", :flip))
	else
		rng = unique(cat(RayPath[1,:], rng, dims = 1))
		plot!(pt, RayPath, vars = (1,2))
	end
	if n == length(θ₀_vec)
		sort!(rng)
		plot!(pt, rng, Bty.z.(rng))
		display(pt)

		rng = RayPath[1,:]
		pt2 = plot(rng, Bty.z.(rng) .- RayPath[2,:])
		display(pt2)
	end
end
end
