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
	dcz
	function Medium(c)
		dcz(z) = ForwardDiff.derivative(c, z)
		return new(c, dcz)
	end
end

"""
TODO:
* Add ray pressure calculation.
* Check initial depth is between bathymetry and altimetry.
* Multiple dispatch on struct contructor methods.
* Report issue to DiffEqBase.
* Generalise inputs.
* Replicate Jensen ray tracing.
"""
function RayTrace(θ₀, Ocn, Src, Bty, Ati, T)
		ξ₀ = cos.(θ₀)/Ocn.c(Src.z)
		ζ₀ = sin.(θ₀)/Ocn.c(Src.z)
		u0 = [Src.r, Src.z, ξ₀, ζ₀]

		function BoundaryReflection!(Ray, Bnd::Boundary)
			θ_bnd = atan(Bnd.dzr(Ray.u[2]))
			θ_inc = asin(Ray.u[4]*Ocn.c(Src.z))
			θ_out = 2θ_bnd - θ_inc
			Ray.u[4] = sin(θ_out)/Ocn.c(Src.z)
			Ray.u[3] = cos(θ_out)/Ocn.c(Src.z)
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
			du[1] = Ocn.c(z)*ξ
			du[2] = Ocn.c(z)*ζ
			du[3] = 0.
			du[4] = -Ocn.dcz(z)/Ocn.c(z)^2
		end

		tspan = (0., T)
		prob = ODEProblem(eikonal!, u0, tspan)

		sol = solve(prob, callback = Bnd_cb)
end

T = 1e4
cMin = 1.5e3
cMax = 1.6e3
cz = [0 0 1; 2.5e5 5e2 1; 1e6 1e3 1.0]\[cMax; cMin; cMax]
cOcn(z) = cMin + ((cz[1]*z^2 + cz[2]*z + cz[3]) - cMin)
r₀ = 0.
z₀ = 250.
zBeg = 1e3
zEnd = 7.5e2
# zBty(r) = zBeg + (zEnd - zBeg)*r/10e3
zBty(r) = 1e3 - 2e2exp(-(r - 5e3)^2/1e6)

Bty = Boundary(zBty)
Ati = Boundary(0.0)
Src = Entity(r₀, z₀)
Ocn = Medium(cOcn)
θ₀_max = acos(cOcn(z₀)/cMax)
θ₀_vec = (-1.5:0.25:1.5).*θ₀_max

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
	end
end
end

##
zBty(r) = 1e3 - 3e2exp(-(r - 5e3)^2/1e6)
dzrBty(r) = ForwardDiff.derivative(zBty, r)
r = range(0, 10e3, length = 100)
plot(r, zBty.(r))
plot(r, dzrBty.(r))