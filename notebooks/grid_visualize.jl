### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ 1c0f68a4-0dd0-11ed-12e9-31130b1d1f6d
begin
	import Pkg; Pkg.activate()
	using Revise
	using GeoPhyInv
	using Statistics
	using Test
	using PlutoUI
	using Luxor
end

# ╔═╡ 2028a25b-7959-4a6f-9a45-94733165c955
md"""
If you are keen on developing finite-difference simulation modules of GeoPhyInv, 
[Luxor graphics](https://juliagraphics.github.io/Luxor.jl/stable/) are available to help 
you easily get started by visualizing 2-D staggered grids. 
Currently, grid visualization in 3D is not available. 
"""

# ╔═╡ 7facc1e5-2946-40d8-a395-68ddcf2f75f5
@init_parallel_stencil(2, true, Float32, 2)

# ╔═╡ 922ef4bb-a6fa-49e7-a608-c942fa294cf1
md"""
Let's get started with plotting 2-D acoustic fields.
"""

# ╔═╡ a4f3f875-30ed-4f0e-8b32-a6243034914e
@draw GeoPhyInv.luxor_mgrid(FdtdAcoustic())


# ╔═╡ 4fe9e85f-b974-435e-a4b8-48fb04bfb5ef
md"""
We can avoid some clutter by isolating fields of interest.
"""

# ╔═╡ 44f19a94-26a5-4050-884a-63af4bfe2452
@draw GeoPhyInv.luxor_mgrid(FdtdAcoustic(), [:vx, :vz])

# ╔═╡ f1af6d1a-3019-4a69-bbe3-30b8cd433d56
md"""
On the other hand, if you like it more clumsier, lets visualize all the staggered grids
for elastic simulation.
"""

# ╔═╡ 9e509b0f-62fa-4260-85d1-708e8ef3c2eb
@draw GeoPhyInv.luxor_mgrid(FdtdElastic())

# ╔═╡ d87e4519-38cb-4c9e-99b5-7849da1db47a
md"""
Lets get rid of those spatial derivatives, and plot only only stress and velocity fields.
"""

# ╔═╡ e0f99e0a-aa32-435a-9948-ab32729612bf
@draw GeoPhyInv.luxor_mgrid(FdtdElastic(), [:tauxx, :tauxz, :tauzz, :vx, :vz])


# ╔═╡ Cell order:
# ╠═2028a25b-7959-4a6f-9a45-94733165c955
# ╠═1c0f68a4-0dd0-11ed-12e9-31130b1d1f6d
# ╠═7facc1e5-2946-40d8-a395-68ddcf2f75f5
# ╠═922ef4bb-a6fa-49e7-a608-c942fa294cf1
# ╠═a4f3f875-30ed-4f0e-8b32-a6243034914e
# ╠═4fe9e85f-b974-435e-a4b8-48fb04bfb5ef
# ╠═44f19a94-26a5-4050-884a-63af4bfe2452
# ╠═f1af6d1a-3019-4a69-bbe3-30b8cd433d56
# ╠═9e509b0f-62fa-4260-85d1-708e8ef3c2eb
# ╠═d87e4519-38cb-4c9e-99b5-7849da1db47a
# ╠═e0f99e0a-aa32-435a-9948-ab32729612bf
