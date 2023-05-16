### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ fcf2e69c-e407-41f6-9113-a7a58a9bbae8
import Pkg; Pkg.activate("GeoPhyInv")

# ╔═╡ fc8eb904-107b-11ed-21a1-b71a02c27438
begin
	using Revise
	using GeoPhyInv
	using Statistics
	using Test
	using PlutoUI
	using Plots; gr()
end

# ╔═╡ af185e45-4203-4c90-8513-fc9cb3672f47
@init_parallel_stencil(2, false, Float32, 2)

# ╔═╡ 15c1460c-29cb-4b79-ad44-0c57c2a48040
md"""
## Custom Construction
Let's begin with 2D. 
To construct an instance of Medium, we need ranges for z and x.
"""

# ╔═╡ 8f1a82a0-71d9-4d7a-a3fe-3c023dbf0e08
mgrid = [range(0.0, stop = 10.0, step = 0.1), range(0.0, stop = 30.0, step = 0.2)];

# ╔═╡ 1fb4497c-8120-4dcb-a05b-0a985dc197f2
md"""
Now we will allocate medium parameters on the grid. To do so, 
we will use `[:vp, :rho]` for an acoustic medium, and add `:vs`
for an elastic medium.
"""

# ╔═╡ 2c15663d-e6dd-438e-ae78-8fc6c6023b8c
medium = Medium(mgrid, [:vp, :rho, :vs])

# ╔═╡ 33051873-3d33-4065-99a1-323aff13e140
md"""
As you can see, the medium just got initiated with all zeros. The first step is to decide on the bounds for the medium parameters. They will be used modeling or inversion. We shall now `update!`
the instance.
"""

# ╔═╡ 0ee6330a-da38-40f9-b290-02b07ba94891
vpb = [2100.0, 2200.0]; vsb = [1500, 1700]; rhob = [2100.0, 2300.0];

# ╔═╡ bd61a821-8a69-47d0-b915-262d2fa77589
update!(medium, [:vp, :vs, :rho], [vpb, vsb, rhob]);

# ╔═╡ 166eb01f-e09b-459a-8bbb-24ee62e115f0
medium

# ╔═╡ e704a7f8-93ca-4dd7-a2aa-3574e49e3e0e
md"""
Now it's time to fill in the actual medium parameters. One option is to simply fill the medium with mean values.
"""

# ╔═╡ 15b1c54e-f9b5-4658-bfc1-81298307cafb
fill!(medium)

# ╔═╡ b6cdb950-8935-43d2-a81c-8b03ec03c5c1
md"""
Otherwise, we can manually update parameters of medium.
"""

# ╔═╡ 96b199f4-5b4d-415d-8a09-918951a0a2a7
begin
	medium[:vp] .= 3000.0;
	medium[:vs] .= 2000.0;
	medium
end

# ╔═╡ e5cfee1e-5876-40e0-aab1-19af9f093468
md"""
## Update
The medium can now be "updated" by adding random noise.
"""

# ╔═╡ 99d3ffc1-0bbe-4b09-a20c-41141e3ac5d5
update!(medium, [:vp, :rho], randn_perc = 1.0); plot(medium)

# ╔═╡ a5f5ffe7-136f-4e4e-9ba7-48a6c07f1e07
md"""
## Gallery 
"""

# ╔═╡ Cell order:
# ╠═fcf2e69c-e407-41f6-9113-a7a58a9bbae8
# ╠═fc8eb904-107b-11ed-21a1-b71a02c27438
# ╠═af185e45-4203-4c90-8513-fc9cb3672f47
# ╠═15c1460c-29cb-4b79-ad44-0c57c2a48040
# ╠═8f1a82a0-71d9-4d7a-a3fe-3c023dbf0e08
# ╠═1fb4497c-8120-4dcb-a05b-0a985dc197f2
# ╠═2c15663d-e6dd-438e-ae78-8fc6c6023b8c
# ╠═33051873-3d33-4065-99a1-323aff13e140
# ╠═0ee6330a-da38-40f9-b290-02b07ba94891
# ╠═bd61a821-8a69-47d0-b915-262d2fa77589
# ╠═166eb01f-e09b-459a-8bbb-24ee62e115f0
# ╠═e704a7f8-93ca-4dd7-a2aa-3574e49e3e0e
# ╠═15b1c54e-f9b5-4658-bfc1-81298307cafb
# ╠═b6cdb950-8935-43d2-a81c-8b03ec03c5c1
# ╠═96b199f4-5b4d-415d-8a09-918951a0a2a7
# ╠═e5cfee1e-5876-40e0-aab1-19af9f093468
# ╠═99d3ffc1-0bbe-4b09-a20c-41141e3ac5d5
# ╠═a5f5ffe7-136f-4e4e-9ba7-48a6c07f1e07
