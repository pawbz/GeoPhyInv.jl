### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 0ffd8ee4-1735-4756-befb-c7c10d08eb34
# ╠═╡ show_logs = false
begin
import Pkg
Pkg.add("PlutoLinks")
Pkg.add("PlutoUI")
Pkg.add("PlutoTest")
Pkg.add("Plots")
using PlutoLinks: @revise
using PlutoUI, PlutoTest, Plots
end

# ╔═╡ 86d3f068-a979-42f5-a9e7-138e94c16b38
@bind reload_geophyinv Button("using GeoPhyInv")

# ╔═╡ 7687d367-f9d5-4539-9156-e26d87379f87
# ╠═╡ show_logs = false
begin
	reload_geophyinv
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
	@revise using GeoPhyInv
end

# ╔═╡ d9b71485-8b64-4ad0-a242-dde6300af835
# ╠═╡ show_logs = false
begin
	reload_geophyinv
	@init_parallel_stencil(2, false, Float32, 2)
	using Statistics
end

# ╔═╡ de2f97fa-f64d-4754-bd34-e44dbf13c336
md"""
In order to install `GeoPhyInv` enter these package manager commands in the REPL.
```julia
using Pkg
Pkg.add(PackageSpec(name="GeoPhyInv",url="https://github.com/pawbz/GeoPhyInv.jl.git"))
```
It is necessary to configure GeoPhyInv with a macro `@init_parallel_stencil` before using it. If you need to change this configuration, the julia kernel must be restarted.
```julia
using GeoPhyInv; @init_parallel_stencil(⋯)
```
"""

# ╔═╡ a6b59ba1-ea71-467f-bf26-49a634a1bbd9
md"""
## AGeom
"""

# ╔═╡ 3c549ecd-52dc-47d2-b846-620148089296
@doc AGeom

# ╔═╡ addd0bef-5d52-4e06-9350-f6dade3db428
md"## Examples"

# ╔═╡ 50494fe6-6bdf-4db1-8e61-6c9b39627b1f
md"""
Let's create a 2-D `mgrid` for the experiment.
"""

# ╔═╡ ae4201a4-7bb1-456c-bdec-d55127324118
mgrid2D = [range(0.0, stop = 10.0, step = 0.1), range(0.0, stop = 30.0, step = 0.2)];

# ╔═╡ a7d1f59b-8ca1-4b32-b3cd-eab98a015bc5
md"""
We can simply initialize an acquisition on `mgrid`, where the positions will be randomly chosen. 
"""

# ╔═╡ c63edf8b-64f5-40ed-a3c0-fa07823c52b1
ageom2D = AGeom(mgrid2D, SSrcs(2), Srcs(10), Recs(10))

# ╔═╡ 8292ec17-0f70-4d7b-a1a0-28c9bbb57166
plot(ageom2D, SSrcs()); plot!(ageom2D, Recs(0))

# ╔═╡ 5d045903-cd72-416d-b5d3-65177bc78d81
md"The acquisition geometry above has 2 supersources, 10 sources, and 10 receiver. Similarly, we can redo for a 3-D grid."

# ╔═╡ d66049f2-28d9-4a8f-9c05-b085af64a9b5
begin
	mgrid3D = fill(range(-10, stop=10, step=0.01),3);
	ageom3D = AGeom(mgrid3D, SSrcs(1), Srcs(10), Recs(10));
end

# ╔═╡ 9a798cae-5cb0-4508-ad57-db50e07f2055
plot(ageom3D, SSrcs()); plot!(ageom3D, Recs(0))

# ╔═╡ f446b1cd-b1c6-4da0-907e-44028ca858ba
md"""
In the case of the 2-D grid, we can also use one of the predefined acquisitions.
"""

# ╔═╡ 07cdf525-94a4-495f-ad7d-14707424a8f4
ageom_xwell = AGeom(mgrid2D, :xwell, SSrcs(3), Recs(10))

# ╔═╡ fe64563a-582d-4ff8-be50-4e2541cdd1f4
plot(ageom_xwell, SSrcs()); plot!(ageom_xwell, Recs(0))

# ╔═╡ f5f697fb-a465-4ce3-b394-03025f0a88bd
md"""
Obviously, the sources and receivers are placed inside mgrid, let's test that.
"""

# ╔═╡ f5022af5-8342-45d8-937a-d2372909960b
@test (ageom2D ∈ mgrid2D)

# ╔═╡ 4c7845e5-eb40-49af-959c-e01375345084
@test (ageom3D ∈ mgrid3D)

# ╔═╡ bf7ca660-16ce-4544-94f9-317ff2568279
md"The source and receiver positions can be updated as desired."

# ╔═╡ 1fb65b26-49be-4f22-8c93-3b63fe6dd8f8
begin
	update!(ageom2D[1], Srcs(1), [0, 1], [10, 20]);
	update!(ageom2D[1], Recs(1), [0, 0], [10, 20]);
	update!(ageom2D, SSrcs(), [0, 1], [10, 20]);
end

# ╔═╡ 49e504be-6a44-45ce-92bf-1bae9bd0f25e
md"It is easy to combine supersources. Now ageom2 has 20 supersources."

# ╔═╡ 2a1a059e-cf52-4e6c-a2dc-a700e4419594
ageom2_new = vcat(ageom2D, ageom_xwell);

# ╔═╡ 636eac63-09c3-492e-add0-55552a9fec00
plot(ageom2_new, SSrcs())

# ╔═╡ Cell order:
# ╟─86d3f068-a979-42f5-a9e7-138e94c16b38
# ╟─de2f97fa-f64d-4754-bd34-e44dbf13c336
# ╟─0ffd8ee4-1735-4756-befb-c7c10d08eb34
# ╟─7687d367-f9d5-4539-9156-e26d87379f87
# ╠═d9b71485-8b64-4ad0-a242-dde6300af835
# ╟─a6b59ba1-ea71-467f-bf26-49a634a1bbd9
# ╠═3c549ecd-52dc-47d2-b846-620148089296
# ╟─addd0bef-5d52-4e06-9350-f6dade3db428
# ╟─50494fe6-6bdf-4db1-8e61-6c9b39627b1f
# ╠═ae4201a4-7bb1-456c-bdec-d55127324118
# ╟─a7d1f59b-8ca1-4b32-b3cd-eab98a015bc5
# ╠═c63edf8b-64f5-40ed-a3c0-fa07823c52b1
# ╠═8292ec17-0f70-4d7b-a1a0-28c9bbb57166
# ╟─5d045903-cd72-416d-b5d3-65177bc78d81
# ╠═d66049f2-28d9-4a8f-9c05-b085af64a9b5
# ╠═9a798cae-5cb0-4508-ad57-db50e07f2055
# ╟─f446b1cd-b1c6-4da0-907e-44028ca858ba
# ╠═07cdf525-94a4-495f-ad7d-14707424a8f4
# ╠═fe64563a-582d-4ff8-be50-4e2541cdd1f4
# ╟─f5f697fb-a465-4ce3-b394-03025f0a88bd
# ╠═f5022af5-8342-45d8-937a-d2372909960b
# ╠═4c7845e5-eb40-49af-959c-e01375345084
# ╠═bf7ca660-16ce-4544-94f9-317ff2568279
# ╠═1fb65b26-49be-4f22-8c93-3b63fe6dd8f8
# ╠═49e504be-6a44-45ce-92bf-1bae9bd0f25e
# ╠═2a1a059e-cf52-4e6c-a2dc-a700e4419594
# ╠═636eac63-09c3-492e-add0-55552a9fec00
