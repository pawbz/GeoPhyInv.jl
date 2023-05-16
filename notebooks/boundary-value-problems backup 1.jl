### A Pluto.jl notebook ###
# v0.19.13

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

# ╔═╡ 45656ece-32e9-488f-be20-6546017e1e94
TableOfContents()

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
	GeoPhyInv.@init_parallel_stencil(2, false, Float32, 2)
	using Statistics, LossFunctions, LinearAlgebra
end

# ╔═╡ de2f97fa-f64d-4754-bd34-e44dbf13c336
md"""
In order to install `GeoPhyInv` enter these package manager commands in the REPL.
```julia
using Pkg
Pkg.add(PackageSpec(name="GeoPhyInv",url="https://github.com/pawbz/GeoPhyInv.jl.git"))
```
It is important to configure GeoPhyInv with a macro `@init_parallel_stencil` before using it. If you need to change this configuration, the julia kernel must be restarted.
```julia
using GeoPhyInv; @init_parallel_stencil(⋯)
```
"""

# ╔═╡ 2c5e29cf-67cd-4c74-b7e8-8f16b1828390
pa = SeisForwExpt(:acou_homo2D);

# ╔═╡ b236913b-5320-461a-8583-74eb4140ff27


# ╔═╡ 4a99d557-2f8b-4b28-9eb6-7aab01c32c27
for field in [:p, :vx, :vz]
    println("############ Testing Backprop for source type ", field)

    srcwav = pa[:srcwav][1]
    GeoPhyInv.update!(srcwav[1], [field])

    for src_types in [[1, -1], [2, -2]]
		# save boundary and final states
        pa.c.backprop_flag = :save

        GeoPhyInv.update!(pa, [srcwav], [src_types[1]])

        update!(pa)
        rec1 = deepcopy(pa.c.data[1])
        rec1 = rec1[1].d[1]
        rec1 ./= std(rec1)

        # change source flag and update wavelets in pa
        GeoPhyInv.update!(pa, [srcwav], [src_types[2]])

		# force boundary values and use saved initial state
        pa.c.backprop_flag = :force 

        update!(pa)
        rec2 = deepcopy(pa.c.data[1])

        # time reverse records before computing L2dist
        reverse!(rec2)
        rec2 = rec2[1].d[1]
        rec2 ./= std(rec2)

        # compute L2dist
        @show err = value(L2DistLoss(), rec1, rec2, AggMode.Mean())

        # desired accuracy?
        @test err < 1e-10
    end
end

# ╔═╡ Cell order:
# ╠═45656ece-32e9-488f-be20-6546017e1e94
# ╟─86d3f068-a979-42f5-a9e7-138e94c16b38
# ╟─de2f97fa-f64d-4754-bd34-e44dbf13c336
# ╟─0ffd8ee4-1735-4756-befb-c7c10d08eb34
# ╟─7687d367-f9d5-4539-9156-e26d87379f87
# ╠═d9b71485-8b64-4ad0-a242-dde6300af835
# ╠═2c5e29cf-67cd-4c74-b7e8-8f16b1828390
# ╠═4a99d557-2f8b-4b28-9eb6-7aab01c32c27
