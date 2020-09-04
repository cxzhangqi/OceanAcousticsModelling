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
	function BoundaryReflection(ξᵢ, ζᵢ, θᵦ, c)
		θᵢ = acos(c*ξᵢ)
		assert(θᵢ == asin(c*ζᵢ))
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
zAti = 0.

# ╔═╡ f4c0f360-eea2-11ea-255d-7918dac724a5
zBty = 1e3

# ╔═╡ 1021ed40-ee98-11ea-332e-33df92f557f0
begin
	C = [1 zAti zAti^2
		 1 (zAti + zBty)/2 ((zAti + zBty)/2)^2
		 1 zBty zBty^2]
	c₀, c₁, c₂ = C\[cMax, cMin, cMax]
	c(r, z) = c₀ + c₁*z + c₂*z^2
	c(x) = c(x[1], x[2])
end

# ╔═╡ 1bc16300-eea3-11ea-2e16-fbed419f1b4a
begin
	zTemp = range(0., 1e3, length = 100)
	plot(c.(0, zTemp), zTemp,
		xaxis = ("Celerity (m/s)"),
		yaxis = ("Depth (m)", :flip))
end

# ╔═╡ 3eb65c80-eea3-11ea-2dbe-e1c9c19f7824
begin
	∇c(x) = ForwardDiff.gradient(c, x)
	∇c(r, z) = ∇c([r, z])
	∂c∂r(r, z) = ∇c(r, z)[1]
	∂c∂z(r, z) = ∇c(r, z)[2]
	
	plot(∂c∂z.(0, zTemp), zTemp,
		yaxis = ("Depth (m)", :flip),
		xaxis = ("∂c/∂z (m/s)"))
end

# ╔═╡ d9811810-ee97-11ea-3cc9-ed4044180fc7
function eikonal!(du, u, p, t)
	r = u[1]
	z = u[2]
	ξ = u[3]
	ζ = u[4]
	dr = c(r, z)*ξ
	dz = c(r, z)*ζ
	dξ = -∂c∂r(r, z)/c(r, z)^2
	dζ = -∂c∂z(r, z)/c(r, z)^2
	du[1] = dr
	du[2] = dz
	du[3] = dξ
	du[4] = dζ
end

# ╔═╡ 8daa55b0-eea0-11ea-08d5-e9557ab3a90b
begin
	r₀ = 0
	z₀ = (zAti + zBty)/2
	θ₀ = 1.2acos(c(r₀, z₀)/cMax)
end

# ╔═╡ cfa55b40-eea0-11ea-211b-2d9c41d749d1
begin
	ξ₀ = cos(θ₀)/c(r₀, z₀)
	ζ₀ = sin(θ₀)/c(r₀, z₀)
	u₀ = [r₀, z₀, ξ₀, ζ₀]
end

# ╔═╡ 22d9ae80-eea4-11ea-36fe-25b56eda117b
begin
	Aty_condition(u, t, integrator) = u[2] - zAti
	Aty_affect!(integrator) = integrator.u[4] *= -1
	Aty_cb = ContinuousCallback(Aty_condition, Aty_affect!)
	
	Bty_condition(u, t, integrator) = u[2] - zBty
	Bty_affect!(integrator) = integrator.u[4] *= -1
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

# ╔═╡ 8c885900-eea2-11ea-2afd-19f585a75836
plot(sol, vars = (1,2))

# ╔═╡ Cell order:
# ╠═ad471390-ee8c-11ea-0f90-a771df3ba437
# ╠═93015d70-ee86-11ea-0d16-d30409dcf0fa
# ╠═324f8070-ee8a-11ea-071e-638691f85469
# ╠═5aeff990-ee8c-11ea-3199-e36fd1e94941
# ╠═8320abd0-ee8c-11ea-2ee9-952362b9e575
# ╟─e4f3ff92-eea2-11ea-1b97-07090817d618
# ╟─ef984d70-eea2-11ea-3fb1-89712c51f9cb
# ╟─f2886882-eea2-11ea-148f-cd4108df600c
# ╟─f4c0f360-eea2-11ea-255d-7918dac724a5
# ╠═1021ed40-ee98-11ea-332e-33df92f557f0
# ╠═1bc16300-eea3-11ea-2e16-fbed419f1b4a
# ╠═3eb65c80-eea3-11ea-2dbe-e1c9c19f7824
# ╠═d9811810-ee97-11ea-3cc9-ed4044180fc7
# ╠═8daa55b0-eea0-11ea-08d5-e9557ab3a90b
# ╠═cfa55b40-eea0-11ea-211b-2d9c41d749d1
# ╠═22d9ae80-eea4-11ea-36fe-25b56eda117b
# ╠═c560c420-eea1-11ea-0046-7f8917d78848
# ╠═8c885900-eea2-11ea-2afd-19f585a75836
