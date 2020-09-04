### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# ╔═╡ ad471390-ee8c-11ea-0f90-a771df3ba437
begin
	using Plots
	using ForwardDiff
	using DifferentialEquations
end

# ╔═╡ 93015d70-ee86-11ea-0d16-d30409dcf0fa
begin
	function BoundaryReflection(θᵢ, θᵦ)
		return θᵣ = 2θᵦ - θᵢ
	end
	function BoundaryReflection(ξᵢ, ζᵢ, ∂zBty∂r, c)
		θᵢ = acos(c*ξᵢ)
		θᵦ = atan(∂zBty∂r)
		if θᵢ ≠ asin(c*ζᵢ)
			ErrorException("ξ and ζ not defining arc-length component.")
		end
		θᵣ = BoundaryReflection(θᵢ, θᵦ)
		ξᵣ = cos(θᵣ)/c
		ζᵣ = sin(θᵣ)/c
		return ξᵣ, ζᵣ
	end
end

# ╔═╡ 324f8070-ee8a-11ea-071e-638691f85469
begin
	θᵢ = 0
	θᵦ = π/5
	θᵣ = BoundaryReflection(θᵢ, θᵦ)
	rad2deg(θᵣ)
end

# ╔═╡ 5aeff990-ee8c-11ea-3199-e36fd1e94941
begin
	xᵢ = [-sqrt(2)*cos(θᵢ), 0]
	yᵢ = [sqrt(2)*sin(θᵢ), 0]
	xᵣ = [0, sqrt(2)*cos(θᵣ)]
	yᵣ = [0, sqrt(2)*sin(θᵣ)]
	xᵦ = range(-1, 1, length = 2)
	yᵦ = sin(θᵦ).*xᵦ
end

# ╔═╡ 8320abd0-ee8c-11ea-2ee9-952362b9e575
begin
	plot(aspect_ratio = 1,
		xaxis = ("x", (-1, 1)),
		yaxis = ("y", (-1, 1)))
	plot!(xᵢ, yᵢ)
	plot!(xᵣ, yᵣ)
	plot!(xᵦ, yᵦ)
end

# ╔═╡ e4f3ff92-eea2-11ea-1b97-07090817d618
cMax = 1600.

# ╔═╡ ef984d70-eea2-11ea-3fb1-89712c51f9cb
cMin = 1500.

# ╔═╡ f2886882-eea2-11ea-148f-cd4108df600c
begin
# 	zAti = 0.
	zAtiMin = 100
	zAtiMax = 300
	zAti(r) = zAtiMin + (zAtiMax - zAtiMin)*(sin(r/1e3) + 1.)/2
end

# ╔═╡ f4c0f360-eea2-11ea-255d-7918dac724a5
begin
	rPeak = 4000.
	rMax = 10e3
	zMax = 10e3
	zMin = 8e3
	Aᵣ = (2rPeak/3)^2/log((9zMax - 11zMin)/(10(zMax - zMin)))
	zBty(r) = zMax - (zMax - zMin)*exp(-(r - rPeak)^2/4e5)
end

# ╔═╡ b5d293c0-eeb1-11ea-33c5-0709cebee220
begin
	rTemp = range(0, rMax, length = 234)
	zTemp = range(0., zMax, length = 123)
end

# ╔═╡ dbf7848e-eebe-11ea-1054-617dccc52071
plot(rTemp, zAti.(rTemp), yaxis = (:flip, (0, zMax)))

# ╔═╡ 92e8fd2e-eeb7-11ea-372e-a730c669fe4f
∂zBty∂r(r) = ForwardDiff.derivative(zBty, r)

# ╔═╡ 97f35f00-eeb2-11ea-066e-2bb671bb0e86
plot(rTemp, zBty.(rTemp), yaxis = (:flip, (zTemp[1], zTemp[end])))

# ╔═╡ 1021ed40-ee98-11ea-332e-33df92f557f0
C(r) = [1 zAti(r) zAti(r)^2
	 1 (zAti(r) + zBty(r))/2 ((zAti(r) + zBty(r))/2)^2
	 1 zBty(r) zBty(r)^2]

# ╔═╡ 55c908e0-eeae-11ea-1ed3-9ba223b87361
c_(r) = C(r)\[cMax, cMin, cMax]

# ╔═╡ fcef8550-eead-11ea-11c1-9bc7f701b93a
begin
	c₀(r) = c_(r)[1]
	c₁(r) = c_(r)[2]
	c₂(r) = c_(r)[3]
end

# ╔═╡ 8142dd90-eeac-11ea-1543-1f0507a04677
begin
	c(r, z) = c₀(r) + c₁(r)*z + c₂(r)*z^2
	c(x) = c(x[1], x[2])
end

# ╔═╡ 2fc217a0-eeb2-11ea-1ea3-6de33f348736
begin
	plot()
	plot!(rTemp, c₀.(rTemp))
	plot!(rTemp, c₁.(rTemp))
	plot!(rTemp, c₂.(rTemp))
end

# ╔═╡ 24825bd0-eeac-11ea-18f7-33e12b8812de
c(300., 500.)

# ╔═╡ d6c61f50-eeae-11ea-2e5e-c9d16b7bd6fc
contour(rTemp, zTemp, c, fill = true, yaxis = :flip)

# ╔═╡ 1bc16300-eea3-11ea-2e16-fbed419f1b4a
begin
	rVal = 3000.
	plot(c.(rVal, zTemp), zTemp,
		xaxis = ("Celerity (m/s)"),
		yaxis = ("Depth (m)", :flip, (zAti(rVal), zBty(rVal))))
end

# ╔═╡ 3eb65c80-eea3-11ea-2dbe-e1c9c19f7824
begin
	∇c(x) = ForwardDiff.gradient(c, x)
	∇c(r, z) = ∇c([r, z])
	∂c∂r(r, z) = ∇c(r, z)[1]
	∂c∂z(r, z) = ∇c(r, z)[2]
end

# ╔═╡ ce610ee0-eeb5-11ea-1903-55211bcd508b
begin
	prr = plot(∂c∂z.(rMax/2, zTemp), zTemp,
		yaxis = ("Depth (m)", :flip),
		xaxis = ("∂c/∂z (m/s)"),
		title = "r = $(rMax/2)")
	prz = plot(rTemp,∂c∂r.(rTemp, zMax/2),
		xaxis = ("Range (m)"),
		yaxis = ("∂c/∂r (m/s)"),
		title = "z = $(zMax/2)")
	pzr = plot(rTemp,∂c∂z.(rTemp, zMax/2),
		xaxis = ("Range (m)"),
		yaxis = ("∂c/∂z (m/s)"),
		title = "z = $(zMax/2)")
	pzz = plot(∂c∂z.(rMax/2, zTemp), zTemp,
		yaxis = ("Depth (m)", :flip),
		xaxis = ("∂c/∂z (m/s)"),
		title = "r = $(rMax/2)")
	l = @layout [a b; c d]
	plot(prr, prz, pzz, pzr, layout = l)
end

# ╔═╡ d9811810-ee97-11ea-3cc9-ed4044180fc7
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

# ╔═╡ 8daa55b0-eea0-11ea-08d5-e9557ab3a90b
begin
	r₀ = 0
	z₀ = (zAti + zBty(r₀))/2
	θ₀ = 2acos(c(r₀, z₀)/cMax)
end

# ╔═╡ cfa55b40-eea0-11ea-211b-2d9c41d749d1
begin
	ξ₀ = cos(θ₀)/c(r₀, z₀)
	ζ₀ = sin(θ₀)/c(r₀, z₀)
	τ₀ = 0
	u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀]
end

# ╔═╡ 22d9ae80-eea4-11ea-36fe-25b56eda117b
begin
	Aty_condition(u, t, integrator) = u[2] - zAti
	Aty_affect!(integrator) = integrator.u[4] *= -1
	Aty_cb = ContinuousCallback(Aty_condition, Aty_affect!)
	
	Bty_condition(u, t, integrator) = u[2] - zBty(u[1])
	Bty_affect!(ray) = ray.u[3], ray.u[4] = BoundaryReflection(ray.u[3], ray.u[4], ∂zBty∂r(ray.u[1]), c(ray.u[1], ray.u[2]))
	Bty_cb = ContinuousCallback(Bty_condition, Bty_affect!)
	
	Bnd_cb = CallbackSet(Aty_cb, Bty_cb)
end

# ╔═╡ c560c420-eea1-11ea-0046-7f8917d78848
begin
	T = 1e4
	tspan = (0., T)
	prob = ODEProblem(eikonal!, u₀, tspan)
	sol = solve(prob, callback = Bnd_cb)
end

# ╔═╡ 5f742150-eebc-11ea-2dfe-5b1750068f5c
struct Ray
	sMax
	τ
	r
	z
	ξ
	ζ
	function Ray(sol)
		sMax = sol.t[end]
		τ(s) = sol(s, idxs = 5)
		r(s) = sol(s, idxs = 1)
		z(s) = sol(s, idxs = 2)
		ξ(s) = sol(s, idxs = 3)
		ζ(s) = sol(s, idxs = 4)
		return new(sMax, τ, r, z, ξ, ζ)
	end
end

# ╔═╡ fce5002e-eebc-11ea-1368-4d4177a0017e
ray = Ray(sol)

# ╔═╡ 8c885900-eea2-11ea-2afd-19f585a75836
begin
	plot(sol, vars = (1,2), yaxis = :flip, aspect_ratio = 1)
	plot!(rTemp, zBty.(rTemp), xaxis = (0, rMax))
# 	sTemp = range(0, ray.sMax, length = 100)
# 	plot(ray.r.(sTemp), ray.z.(sTemp), yaxis = :flip)
# 	plot!(ray.r.(sTemp), zBty.(ray.r.(sTemp)))
end

# ╔═╡ 3519da50-eea6-11ea-1925-01e1b341d64b
plot(sol, vars = (0, 5))

# ╔═╡ Cell order:
# ╠═ad471390-ee8c-11ea-0f90-a771df3ba437
# ╠═93015d70-ee86-11ea-0d16-d30409dcf0fa
# ╠═324f8070-ee8a-11ea-071e-638691f85469
# ╠═5aeff990-ee8c-11ea-3199-e36fd1e94941
# ╠═8320abd0-ee8c-11ea-2ee9-952362b9e575
# ╟─e4f3ff92-eea2-11ea-1b97-07090817d618
# ╟─ef984d70-eea2-11ea-3fb1-89712c51f9cb
# ╠═b5d293c0-eeb1-11ea-33c5-0709cebee220
# ╠═f2886882-eea2-11ea-148f-cd4108df600c
# ╠═dbf7848e-eebe-11ea-1054-617dccc52071
# ╠═f4c0f360-eea2-11ea-255d-7918dac724a5
# ╠═92e8fd2e-eeb7-11ea-372e-a730c669fe4f
# ╠═97f35f00-eeb2-11ea-066e-2bb671bb0e86
# ╠═1021ed40-ee98-11ea-332e-33df92f557f0
# ╠═55c908e0-eeae-11ea-1ed3-9ba223b87361
# ╠═fcef8550-eead-11ea-11c1-9bc7f701b93a
# ╠═8142dd90-eeac-11ea-1543-1f0507a04677
# ╠═2fc217a0-eeb2-11ea-1ea3-6de33f348736
# ╠═24825bd0-eeac-11ea-18f7-33e12b8812de
# ╠═d6c61f50-eeae-11ea-2e5e-c9d16b7bd6fc
# ╠═1bc16300-eea3-11ea-2e16-fbed419f1b4a
# ╠═3eb65c80-eea3-11ea-2dbe-e1c9c19f7824
# ╠═ce610ee0-eeb5-11ea-1903-55211bcd508b
# ╠═d9811810-ee97-11ea-3cc9-ed4044180fc7
# ╠═8daa55b0-eea0-11ea-08d5-e9557ab3a90b
# ╠═cfa55b40-eea0-11ea-211b-2d9c41d749d1
# ╠═22d9ae80-eea4-11ea-36fe-25b56eda117b
# ╠═c560c420-eea1-11ea-0046-7f8917d78848
# ╠═5f742150-eebc-11ea-2dfe-5b1750068f5c
# ╠═fce5002e-eebc-11ea-1368-4d4177a0017e
# ╠═8c885900-eea2-11ea-2afd-19f585a75836
# ╠═3519da50-eea6-11ea-1925-01e1b341d64b
