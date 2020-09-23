function boundary_reflection_vectors(t_inc::Vector, t_bnd::Vector)
	# works for generic boundary
	n_bnd = [-t_bnd[2], t_bnd[1]]
# 	t_rfl = t_inc - 2(t_inc ⋅ n_bnd)*n_bnd
	t_rfl = t_inc - 2LinearAlgebra.dot(t_inc, n_bnd)*n_bnd

	MyAngle(tng) = atand(tng[2]/tng[1])
	θ_inc = MyAngle(t_inc)
	θ_bnd = MyAngle(t_bnd)
	θ_rfl = MyAngle(t_rfl)
	println(θ_inc)
	println(θ_bnd)
	println(θ_rfl)

	return t_rfl
end

function boundary_reflection_angles(t_inc::Vector, t_bnd::Vector)
	# works for parabolic boundary
	MyAngle(tng) = atan(tng[2]/tng[1])
	θ_inc = MyAngle(t_inc)
	θ_bnd = MyAngle(t_bnd)

	θ_inc_flat = θ_inc - θ_bnd
	θ_rfl_flat = -θ_inc_flat
	θ_rfl = θ_rfl_flat + θ_bnd
	return [cos(θ_rfl), sin(θ_rfl)]
end

θ_inc = π/4
θ_bnd = π/8
boundary_reflection_vectors()