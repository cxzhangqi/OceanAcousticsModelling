## Ray Tracing
using DifferentialEquations

cMin = 1500
cMax = 1600
α = [0 0 1; 2.5e5 5e2 1; 1e6 1e3 1]\[cMax; cMin; cMax]
c(r, z) = α[1]*z^2 + α[2]*z + α[3]
# c(r, z) = 0. < z < 1e3 ? c_(r, z) : -c_(r, z)
∂cz(r, z) = 2α[1]*z + α[2]
∂cr(r, z) = 0.

# Depth = range(0., 1e3, length = 100)
# Speed = c.(0, Depth)

# plot(Speed, Depth)

r₀ = 0.
z₀ = 750.
# θ₀ = atan(10.0/9.1)
θ₀_max = acos(c(r₀, z₀)/cMax)
θ₀_vec = [-θ₀_max; 0; θ₀_max]

tspan = (0., 1e4)
let pt
for n = 1:length(θ₀_vec)
	θ₀ = θ₀_vec[n]
	ξ₀ = cos.(θ₀)/c(r₀, z₀)
	ζ₀ = sin.(θ₀)/c(r₀, z₀)
	u0 = [r₀, z₀, ξ₀, ζ₀]

	function eikonal!(du, u, p, t)
		r = u[1]
		z = u[2]
		ξ = u[3]
		ζ = u[4]
		du[1] = c(r, z)*ξ
		du[2] = c(r, z)*ζ
		du[3] = -∂cr(r, z)/c(r, z)^2
		du[4] = -∂cz(r, z)/c(r, z)^2
	end

	prob = ODEProblem(eikonal!, u0, tspan)

	sol = solve(prob)
	if n == 1
		pt = plot(sol, vars = (1,2))
	else
		plot!(pt, sol, vars = (1,2))
		if n == length(θ₀_vec)
			display(pt)
		end
	end
end
end

## Add Boundaries to Ray Trace
using DifferentialEquations

cMin = 1500
cMax = 1600
α = [0 0 1; 2.5e5 5e2 1; 1e6 1e3 1]\[cMax; cMin; cMax]
c(r, z) = α[1]*z^2 + α[2]*z + α[3]
# c(r, z) = 0. < z < 1e3 ? c_(r, z) : -c_(r, z)
∂cz(r, z) = 2α[1]*z + α[2]
∂cr(r, z) = 0.

# Depth = range(0., 1e3, length = 100)
# Speed = c.(0, Depth)

# plot(Speed, Depth)

r₀ = 0.
z₀ = 250.
# θ₀ = atan(10.0/9.1)
θ₀_max = acos(c(r₀, z₀)/cMax)
# θ₀_vec = [-θ₀_max; 0; θ₀_max]
θ₀_vec = [-1.5; -1.0; 0.0; 1.0; 1.5].*θ₀_max

# Boundaries
# * Altimetry
Aty_condition(u, t, integrator) = u[2]
Aty_affect!(integrator) = integrator.u[4] *= -1
Aty_cb = ContinuousCallback(Aty_condition, Aty_affect!)

# * Bathymetry
z_bty = 1e3
Bty_condition(u, t, integrator) = u[2] - z_bty
Bty_affect!(integrator) = integrator.u[4] *= -1
Bty_cb = ContinuousCallback(Bty_condition, Bty_affect!)

Bnd_cb = CallbackSet(Aty_cb, Bty_cb)

tspan = (0., 1e4)
let pt
for n = 1:length(θ₀_vec)
	θ₀ = θ₀_vec[n]
	ξ₀ = cos.(θ₀)/c(r₀, z₀)
	ζ₀ = sin.(θ₀)/c(r₀, z₀)
	u0 = [r₀, z₀, ξ₀, ζ₀]

	function eikonal!(du, u, p, t)
		r = u[1]
		z = u[2]
		ξ = u[3]
		ζ = u[4]
		du[1] = c(r, z)*ξ
		du[2] = c(r, z)*ζ
		du[3] = -∂cr(r, z)/c(r, z)^2
		du[4] = -∂cz(r, z)/c(r, z)^2
	end

	prob = ODEProblem(eikonal!, u0, tspan)

	sol = solve(prob, callback = Bnd_cb)
	if n == 1
		pt = plot(sol, vars = (1,2))
	else
		plot!(pt, sol, vars = (1,2))
		if n == length(θ₀_vec)
			display(pt)
		end
	end
end
end

## Encapsulate
using DifferentialEquations, Plots

"""
	eikonal!(du, u, p, t)

Eikonal equation for ray tracing.
"""
function eikonal!(du, u, p, t)
	r = u[1]
	z = u[2]
	ξ = u[3]
	ζ = u[4]
	du[1] = c(r, z)*ξ
	du[2] = c(r, z)*ζ
	du[3] = -∂cr(r, z)/c(r, z)^2
	du[4] = -∂cz(r, z)/c(r, z)^2
end

"""
TODO:
* Check initial depth is between bty and aty.
* Implement automatic differentiation for sound speed profile
"""
function RayTrace(θ₀, c, Src, Bty, Aty)
		ξ₀ = cos.(θ₀)/c(Src.r, Src.z)
		ζ₀ = sin.(θ₀)/c(Src.r, Src.z)
		u0 = [Src.r, Src.z, ξ₀, ζ₀]

		# Boundaries
		# * Altimetry
		Aty_condition(u, t, integrator) = u[2] - Aty.z
		Aty_affect!(integrator) = integrator.u[4] *= -1
		Aty_cb = ContinuousCallback(Aty_condition, Aty_affect!)
		
		# * Bathymetry
		Bty_condition(u, t, integrator) = u[2] - Bty.z
		Bty_affect!(integrator) = integrator.u[4] *= -1
		Bty_cb = ContinuousCallback(Bty_condition, Bty_affect!)
		
		Bnd_cb = CallbackSet(Aty_cb, Bty_cb)

		tspan = (0., 1e4)
		prob = ODEProblem(eikonal!, u0, tspan)

		sol = solve(prob, callback = Bnd_cb)
end

struct Boundary
	z
end

struct Entity
	r
	z
end

struct Medium
	c
	∂cz
	∂cr
end

cMin = 1.5e3
cMax = 1.6e3
cz = [0 0 1; 2.5e5 5e2 1; 1e6 1e3 1.0]\[cMax; cMin; cMax]
cr = log(2)/8e3
c(r, z) = cMin + ((cz[1]*z^2 + cz[2]*z + cz[3]) - cMin)*exp(-cr*r)
∂cz(r, z) = exp(-cr*r)*(2cz[1]*z + cz[2])
∂cr(r, z) = -cr*exp(-cr*r)*((cz[1]*z^2 + cz[2]*z + cz[3]) - cMin)
r₀ = 0.
z₀ = 250.

Bty = Boundary(1e3)
Aty = Boundary(0.0)
Src = Entity(r₀, z₀)
Ocn = Medium(c, ∂cz, ∂cr)
θ₀_max = acos(c(r₀, z₀)/cMax)
θ₀_vec = [-1.5; -1; -0.5; 0.0; 0.5; 1.0; 1.5].*θ₀_max

let pt
for n = 1:length(θ₀_vec)
	θ₀ = θ₀_vec[n]
	RayPath = RayTrace(θ₀, c, Src, Bty, Aty)
	if n == 1
		pt = plot(RayPath,
			vars = (1,2), 
			xaxis = ("Range (m)"),
			yaxis = ("Depth (m)", :flip))
	else
		plot!(pt, RayPath, vars = (1,2))
		if n == length(θ₀_vec)
			display(pt)
		end
	end
end
end

## Range-independent (for now)
using DifferentialEquations, Plots, ForwardDiff

"""
	eikonal!(du, u, p, t)

Eikonal equation for ray tracing.
"""
function eikonal!(du, u, p, t)
	r = u[1]
	z = u[2]
	ξ = u[3]
	ζ = u[4]
	du[1] = c(z)*ξ
	du[2] = c(z)*ζ
	du[3] = 0.
	du[4] = -dcz(z)/c(z)^2
end

"""
TODO:
* Check initial depth is between bty and aty.
* Implement automatic differentiation for sound speed profile
"""
function RayTrace(θ₀, c, Src, Bty, Aty)
		ξ₀ = cos.(θ₀)/c(Src.z)
		ζ₀ = sin.(θ₀)/c(Src.z)
		u0 = [Src.r, Src.z, ξ₀, ζ₀]

		# Boundaries
		# * Altimetry
		Aty_condition(u, t, integrator) = u[2] - Aty.z
		Aty_affect!(integrator) = integrator.u[4] *= -1
		Aty_cb = ContinuousCallback(Aty_condition, Aty_affect!)
		
		# * Bathymetry
		Bty_condition(u, t, integrator) = u[2] - Bty.z
		Bty_affect!(integrator) = integrator.u[4] *= -1
		Bty_cb = ContinuousCallback(Bty_condition, Bty_affect!)
		
		Bnd_cb = CallbackSet(Aty_cb, Bty_cb)

		tspan = (0., 1e4)
		prob = ODEProblem(eikonal!, u0, tspan)

		sol = solve(prob, callback = Bnd_cb)
end

struct Boundary
	z
end

struct Entity
	r
	z
end

struct Medium
	c
	dcz
end

cMin = 1.5e3
cMax = 1.6e3
cz = [0 0 1; 2.5e5 5e2 1; 1e6 1e3 1.0]\[cMax; cMin; cMax]
c(z) = cMin + ((cz[1]*z^2 + cz[2]*z + cz[3]) - cMin)
# dcz(z) = 2cz[1]*z + cz[2]
dcz(z) = ForwardDiff.derivative(c, z)
r₀ = 0.
z₀ = 250.

Bty = Boundary(1e3)
Aty = Boundary(0.0)
Src = Entity(r₀, z₀)
Ocn = Medium(c, dcz)
θ₀_max = acos(c(z₀)/cMax)
# θ₀_vec = [-1.5; -1; -0.5; 0.0; 0.5; 1.0; 1.5].*θ₀_max
θ₀_vec = (-1.5:0.25:1.5).*θ₀_max

let pt
for n = 1:length(θ₀_vec)
	θ₀ = θ₀_vec[n]
	RayPath = RayTrace(θ₀, c, Src, Bty, Aty)
	if n == 1
		pt = plot(RayPath,
			vars = (1,2), 
			xaxis = ("Range (m)"),
			yaxis = ("Depth (m)", :flip))
	else
		plot!(pt, RayPath, vars = (1,2))
		if n == length(θ₀_vec)
			display(pt)
		end
	end
end
end

##
using ForwardDiff

struct Derivatives
	f
	∂f_∂x
	∂f_∂y
	function Derivatives(f)
		f_(z) = f(z[1], z[2])
		∇f(z) = ForwardDiff.gradient(f_, z)
		∇f(x, y) = ∇f([x, y])
		∂f_∂x(x, y) = ∇f(x, y)[1]
		∂f_∂y(x, y) = ∇f(x, y)[2]
		return new(f, ∂f_∂x, ∂f_∂x)
	end
end

f(x,y) = x^2/sqrt(y)
df = Derivatives(f)
df.∂f_∂x(1, 2)

##
z₁ = 1e3
cMin = 1500
cMax = 1600
c(z) = cMin + (cMax - cMin)*abs(z - z₁/2)/(z₁/2)
z = range(0, z₁, length = 101)
plot(z, c.(z))