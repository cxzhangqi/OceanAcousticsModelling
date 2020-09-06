### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# ╔═╡ d2b2ab10-f03c-11ea-2afa-37ea84854c57
begin
	using Plots; plotlyjs()
	using ForwardDiff
	using DifferentialEquations
	using LinearAlgebra
	using ColorSchemes
end

# ╔═╡ 302d7c1e-f03d-11ea-077e-81bad303fe2c
# 1. Inputs
begin
	# Altimetry
	zAtiMin = -10
	zAtiMax = 50
	zAti(r) = zAtiMin + (zAtiMax - zAtiMin)*(sin(r/1e3) + 1.)/2
# 	zAti(r) = 0.
	
	# Bathymetry
	rPeak = 4e3
	rMax = 10e3
	zMax = 1e3
	zMin = 7e2
	Aᵣ = (2rPeak/3)^2/log((9zMax - 11zMin)/(10(zMax - zMin)))
	zBty(r) = zMax - (zMax - zMin)*exp(-(r - rPeak)^2/1e5)
	
	# Ocean
	cMin = 1500
	cMax = 1600
	C(r) = [1 zAti(r) zAti(r)^2
		1 (zAti(r) + zBty(r))/2 ((zAti(r) + zBty(r))/2)^2
		1 zBty(r) zBty(r)^2]
	c_(r) = C(r)\[cMax, cMin, cMax]
	c₀(r) = c_(r)[1]
	c₁(r) = c_(r)[2]
	c₂(r) = c_(r)[3]
	c(r, z) = c₀(r) + c₁(r)*z + c₂(r)*z^2
# 	c(r, z) = cMin + (cMax - cMin)*(z - zAti(r))/(zBty(r) - zAti(r))
	
	# Initial Conditions
	r₀ = 0
	z₀ = (zBty(r₀) + zAti(r₀))/2
# 	θ₀ = 2acos(c(r₀, z₀)/cMax)
	θ₀s = acos(c(r₀, z₀)/cMax).*(-1.5:0.5:1.5)
	
	# Other
	T = 1e4 # Figure out how to replace with condition
end

# ╔═╡ fc809290-f03c-11ea-1697-b15add88e425
begin
	function helmholtz(θ₀, r₀, z₀, c, zAti, zBty, T)
		c_(x) = c(x[1], x[2])
		∇c_(x) = ForwardDiff.gradient(c_, x)
		∇c(r, z) = ∇c_([r, z])
		∂c∂r(r, z) = ∇c(r, z)[1]
		∂c∂z(r, z) = ∇c(r, z)[2]

		function eikonal!(du, u, p, t)
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
		tspan = (0., T)
		prob = ODEProblem(eikonal!, u₀, tspan)
		@time RaySol = solve(prob, callback = cb_bnd)
		return RaySol
	end
	function boundary_reflection(t_inc::Vector, t_bnd::Vector)
		n_bnd = [-t_bnd[2], t_bnd[1]]
		t_rfl = t_inc - 2(t_inc ⋅ n_bnd)*n_bnd
		return t_rfl
	end
end

# ╔═╡ d65a7600-f03f-11ea-069e-f981736436e9
# 2. Solve
@time RaySols = helmholtz.(θ₀s, r₀, z₀, c, zAti, zBty, T)


# ╔═╡ 8f0bb750-f049-11ea-31c5-2baf59172716
#. 3. Plot: Inputs
begin
	r = range(0, rMax, length = 1000)
	z = range(min(0, zAtiMin), zMax, length = 100)
end

# ╔═╡ c1426020-f03f-11ea-261c-319ba22ec9ba
# 3. Plot: Slice
begin
	plot(legend = :outerbottomright,
		xaxis = "Range (m)",
		yaxis = ("Depth", :flip)
	)
	plot!(r, zAti, label = "Altimetry")
	plot!.(RaySols, vars = (1, 2), label = false)
	plot!(r, zBty, label = "Bathymetry")
end

# ╔═╡ 4b40a650-f041-11ea-08c0-ed27ea70a1d5
begin
	# 3. Plot: Celerity
	contour(r, z, c, 
		fill = true,
		title = "Sound Speed Field",
		yaxis = ("Depth", :flip),
		xaxis = "Range (m)",
		seriescolor = cgrad(colorschemes[:jet]),
		clim = (1500, 1600))
	plot!(r, zBty)
end

# ╔═╡ 75e7f480-f05a-11ea-2403-3983dce22a56
# TODO: Implement Sound Speed Curves over discrete ranges

# ╔═╡ Cell order:
# ╠═d2b2ab10-f03c-11ea-2afa-37ea84854c57
# ╠═fc809290-f03c-11ea-1697-b15add88e425
# ╠═302d7c1e-f03d-11ea-077e-81bad303fe2c
# ╠═d65a7600-f03f-11ea-069e-f981736436e9
# ╠═8f0bb750-f049-11ea-31c5-2baf59172716
# ╠═c1426020-f03f-11ea-261c-319ba22ec9ba
# ╠═4b40a650-f041-11ea-08c0-ed27ea70a1d5
# ╠═75e7f480-f05a-11ea-2403-3983dce22a56
