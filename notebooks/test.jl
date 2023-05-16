### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ 92194733-d993-4177-8e7e-a545acfd51f8
# ╠═╡ show_logs = false
begin
import Pkg
Pkg.add("PlutoLinks")
Pkg.add("PlutoUI")
Pkg.add("PlutoTest")
Pkg.add("Plots")
	Pkg.add("CUDA")
using PlutoLinks: @revise
using PlutoUI, PlutoTest, Plots, CUDA
end

# ╔═╡ dc1a4656-ee20-11ed-1af5-1b7420730607
begin
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	@revise using GeoPhyInv
end

# ╔═╡ 180d8c39-8a5e-433a-8132-3171d03c29ff
using Interpolations

# ╔═╡ dc3fb7d8-1df8-462b-a0a8-c67c938e2bae
using LinearAlgebra

# ╔═╡ 7d4d6b8c-f03a-46a3-917b-869c8ea9017f
begin
	xs = 1:0.2:5
	ys = 1:0.5:5
	x1 = randn(size(xs))
	# Create linear interpolation object without extrapolation
	
	interp_linear1 = interpolate((xs,), x1, Gridded(Linear()))#extrapolation_bc=Line())
	# interp_linear(3) # exactly log(3)
	y1  = interp_linear1.(ys) # approximately log(3.1)

	y2 = randn(size(ys))
	interp_linear2 = interpolate((ys,), y2, Gridded(Linear()))#extrapolation_bc=Line())
	x2 = interp_linear2.(xs)
	
	# interp_linear(0.9) # outside grid: error
end

# ╔═╡ 2bcd9290-e973-40db-8c31-685c8f210262
gradient(interp_linear1, 1)

# ╔═╡ d2f67504-6fd3-46ed-b6a7-ea25f95fadf9
dot(y2, y1), dot(x1, x2)

# ╔═╡ Cell order:
# ╠═92194733-d993-4177-8e7e-a545acfd51f8
# ╠═dc1a4656-ee20-11ed-1af5-1b7420730607
# ╠═180d8c39-8a5e-433a-8132-3171d03c29ff
# ╠═7d4d6b8c-f03a-46a3-917b-869c8ea9017f
# ╠═2bcd9290-e973-40db-8c31-685c8f210262
# ╠═dc3fb7d8-1df8-462b-a0a8-c67c938e2bae
# ╠═d2f67504-6fd3-46ed-b6a7-ea25f95fadf9
