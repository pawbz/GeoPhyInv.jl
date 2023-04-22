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

# ╔═╡ 4d0a66b2-3330-42d3-b323-4df139fbae53
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

# ╔═╡ 17a7cd01-f98e-4cb8-a93b-f80c713d5423
using LinearMaps

# ╔═╡ f928f64d-6aeb-4081-848e-57e02a8f1030
using DistributedArrays

# ╔═╡ 88854744-d93d-44aa-ab9d-fc26e66572f1
TableOfContents()

# ╔═╡ e86a58b8-23ba-4720-8239-cd5f144649cd
@bind reload_geophyinv Button("using GeoPhyInv")

# ╔═╡ 970c1837-e646-4a3f-9946-9b494af2c21c
begin
    reload_geophyinv
    Pkg.activate(Base.current_project())
    Pkg.instantiate()
    @revise using GeoPhyInv
end

# ╔═╡ a34fc72e-8f98-4263-b48e-6498fa69f211
begin
    reload_geophyinv
	# if CUDA.functional()
			# GeoPhyInv.@init_parallel_stencil(2, true, Float32, 2)
# else
	GeoPhyInv.@init_parallel_stencil(2, false, Float32, 2)

# end
   
	using Statistics
end

# ╔═╡ 8ae50e30-2dc2-490e-953d-8589302c5d5c
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

# ╔═╡ 78769d02-829e-407a-a211-e0d02928e9a8
medium=Medium(:pizza, 10)

# ╔═╡ aafae9ef-f42c-4a93-a2b5-3de6c59ce8bf
plot(medium)

# ╔═╡ d03a303f-33ef-4d85-8162-092478d512f3
medium0=Medium(:acou_homo2D, 10)

# ╔═╡ 0f2d0963-b91a-4b1c-92d5-31f077e28516
ageom = AGeom(medium.mgrid, :surf, SSrcs(1), Recs(100));

# ╔═╡ d72f7155-b766-46ad-a3af-8390b3cfa1f0
plot(medium0, ageom)

# ╔═╡ 120e51bc-cfd4-4a2a-bcd1-fe942ef224d7
begin
	wav, tgrid = ricker(medium, 10, 2)
	srcwav = SrcWav(tgrid, ageom, [:p]); 
	update!(srcwav, [:p], wav)
end

# ╔═╡ a18fbf29-803c-4f0d-9060-fa1201990b16
plot(srcwav[1])

# ╔═╡ 81a68925-ef50-4eda-aacc-3c5acf519dc7


# ╔═╡ cc8fddd6-fe22-45d7-8abb-e9c66f7c449c
pa_acoustic = SeisForwExpt(
    FdtdAcoustic{Born}(),
    medium = medium0,
    ageom = ageom,
	tgrid = tgrid,
    srcwav = srcwav,
	rfields = [:p],
   );

# ╔═╡ c54ae16f-7c39-4899-99e6-6d82c8929bde
GeoPhyInv.get_update_parameters(pa_acoustic)
# dd::Union{GeoPhyInv.PFdtd{FdtdAcoustic{Born}}, GeoPhyInv.PFdtd{FdtdElastic{Born}}}

# ╔═╡ 0836a1ac-ed96-48d6-8a90-73472e5d72a9
typeof(pa_acoustic)

# ╔═╡ ffe9c80a-0230-44e9-8093-8d80da098050
pa_acoustic.c.mod |> names

# ╔═╡ d4d425ea-a1f7-42d1-8c22-562cb7260f98
begin
	update!(pa_acoustic, medium0, medium);
	update!(pa_acoustic, pa_acoustic.c.srcwav, [2, 2], verbose=false)
	# update!(pa_acoustic, medium);
	update!(pa_acoustic)
	# , activepw=[1,2], src_flags=[true, false], rec_flags=[false, true]);
end

# ╔═╡ ab2903f4-b077-4cca-8f5c-ad13ae68b7b0
heatmap(Array(pa_acoustic.c.δmod[:K]))

# ╔═╡ 75825692-7105-4382-b099-7736eea44af1
heatmap(pa_acoustic.c.gradients[:K])

# ╔═╡ d23a59fe-45cb-4898-9e63-d4a3c61b0992
heatmap(Array(pa_acoustic.c.δmod[:rhoI]))

# ╔═╡ 6e59506f-9cc7-4db4-9c6c-8b4a38499074
F = LinearMap(pa_acoustic)

# ╔═╡ 8a43f760-74f3-4a08-b227-bcfbf2894972
# pa_acoustic.c.δmod[:rhoI][250,250] = 1.

# ╔═╡ ac74b62c-b6f6-4a27-83d3-e03e81e843d0
# x=Iterators.flatten(pa_acoustic.c.δmod)|> collect

# ╔═╡ 1f142c02-659e-43d4-a9c1-2e9aa3e2128e
typeof(pa_acoustic.c.δmod)

# ╔═╡ 2b30489a-d03b-4c3b-a75b-3a764156032e
begin
	x= zeros(size(F, 2)); 
	x[div(length(x), 4)]=10000; 
	y = F * x;
	xx = F' * y
end

# ╔═╡ c9a6db62-ee71-4f62-94d0-de43f66e0342
heatmap(reshape(x[1:243*243], 243, 243))

# ╔═╡ 19d90bda-f51e-4c99-8abe-ea5cbe0b39e8
heatmap(reshape(xx[1:243*243], 243, 243))

# ╔═╡ fd3c3a48-d722-4da4-b558-0599af0d0b1c
heatmap(reshape(y, length(tgrid), :), margin=2Plots.cm)

# ╔═╡ 77e288e1-3c98-4268-a243-99c76441520e
plot(pa_acoustic[:data, 1], 1)

# ╔═╡ ee3e5650-553e-4206-b5ab-8c56d7024d1b
iszero(localpart(pa_acoustic.p)[1].ss[1].gradients[:K])

# ╔═╡ 0af23a36-b664-4481-8330-d01eed276025
heatmap((localpart(pa_acoustic.p)[2].w1[:t][:p]))

# ╔═╡ 01d02b35-52aa-48ab-91ab-9ecc6c080043
vec(pa_acoustic.c.mod)

# ╔═╡ e890051d-9c60-4cb7-8f3a-e35fed263bfd
typeof(pa_acoustic[:data, 1])

# ╔═╡ 249bd1ea-e846-4ed5-b52c-bd4f9dd257fe
length(pa_acoustic[:data, 1])

# ╔═╡ 787d4198-cda8-44f0-b629-ed4fd491d9b6
vec(pa_acoustic[:data, 1])

# ╔═╡ 162c3462-e734-4ab6-bc7c-b4fc1f7d763f
heatmap(pa_acoustic.c.data[2][1][:p])

# ╔═╡ c14e3875-bc4a-46bf-8476-371628fccb51
pa_acoustic.c.data[2][1].d

# ╔═╡ e79340da-84ac-4cb3-a90a-49058f541597
pa_acoustic.c.attrib_mod.mode

# ╔═╡ f9e49031-b95d-4578-9360-4fa645b66459
broadcast(ageom) do a
	AGeomss(a.r, a.r, a.nr, a.nr)
end

# ╔═╡ 2ddf7432-4774-404d-bfea-5b90dcf7f7a6
srcwav[1].d |> names

# ╔═╡ 62054879-adc0-4dc2-8bbf-9ab7019a8cd6
Iterators.flatten(pa_acoustic.c.data[2][1].d)

# ╔═╡ f3046085-e49e-4e50-b93e-8a3969518852
pa_acoustic.c.δmod[:rhoI]

# ╔═╡ 5a0e79d2-93e9-4b8e-8ef5-9045a0bf01ad
medium[:K]

# ╔═╡ Cell order:
# ╟─88854744-d93d-44aa-ab9d-fc26e66572f1
# ╟─e86a58b8-23ba-4720-8239-cd5f144649cd
# ╟─8ae50e30-2dc2-490e-953d-8589302c5d5c
# ╠═4d0a66b2-3330-42d3-b323-4df139fbae53
# ╟─970c1837-e646-4a3f-9946-9b494af2c21c
# ╠═a34fc72e-8f98-4263-b48e-6498fa69f211
# ╠═78769d02-829e-407a-a211-e0d02928e9a8
# ╠═aafae9ef-f42c-4a93-a2b5-3de6c59ce8bf
# ╠═d03a303f-33ef-4d85-8162-092478d512f3
# ╠═d72f7155-b766-46ad-a3af-8390b3cfa1f0
# ╠═120e51bc-cfd4-4a2a-bcd1-fe942ef224d7
# ╠═0f2d0963-b91a-4b1c-92d5-31f077e28516
# ╠═c54ae16f-7c39-4899-99e6-6d82c8929bde
# ╠═a18fbf29-803c-4f0d-9060-fa1201990b16
# ╠═81a68925-ef50-4eda-aacc-3c5acf519dc7
# ╠═cc8fddd6-fe22-45d7-8abb-e9c66f7c449c
# ╠═0836a1ac-ed96-48d6-8a90-73472e5d72a9
# ╠═ffe9c80a-0230-44e9-8093-8d80da098050
# ╠═d4d425ea-a1f7-42d1-8c22-562cb7260f98
# ╠═ab2903f4-b077-4cca-8f5c-ad13ae68b7b0
# ╠═75825692-7105-4382-b099-7736eea44af1
# ╠═d23a59fe-45cb-4898-9e63-d4a3c61b0992
# ╠═17a7cd01-f98e-4cb8-a93b-f80c713d5423
# ╠═6e59506f-9cc7-4db4-9c6c-8b4a38499074
# ╠═8a43f760-74f3-4a08-b227-bcfbf2894972
# ╠═ac74b62c-b6f6-4a27-83d3-e03e81e843d0
# ╠═1f142c02-659e-43d4-a9c1-2e9aa3e2128e
# ╠═2b30489a-d03b-4c3b-a75b-3a764156032e
# ╠═c9a6db62-ee71-4f62-94d0-de43f66e0342
# ╠═19d90bda-f51e-4c99-8abe-ea5cbe0b39e8
# ╠═fd3c3a48-d722-4da4-b558-0599af0d0b1c
# ╠═77e288e1-3c98-4268-a243-99c76441520e
# ╠═f928f64d-6aeb-4081-848e-57e02a8f1030
# ╠═ee3e5650-553e-4206-b5ab-8c56d7024d1b
# ╠═0af23a36-b664-4481-8330-d01eed276025
# ╠═01d02b35-52aa-48ab-91ab-9ecc6c080043
# ╠═e890051d-9c60-4cb7-8f3a-e35fed263bfd
# ╠═249bd1ea-e846-4ed5-b52c-bd4f9dd257fe
# ╠═787d4198-cda8-44f0-b629-ed4fd491d9b6
# ╠═162c3462-e734-4ab6-bc7c-b4fc1f7d763f
# ╠═c14e3875-bc4a-46bf-8476-371628fccb51
# ╠═e79340da-84ac-4cb3-a90a-49058f541597
# ╠═f9e49031-b95d-4578-9360-4fa645b66459
# ╠═2ddf7432-4774-404d-bfea-5b90dcf7f7a6
# ╠═62054879-adc0-4dc2-8bbf-9ab7019a8cd6
# ╠═f3046085-e49e-4e50-b93e-8a3969518852
# ╠═5a0e79d2-93e9-4b8e-8ef5-9045a0bf01ad
