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

# ╔═╡ e2e1d696-ad9d-4e87-9e9a-11a3f68d05de
# ╠═╡ show_logs = false
Pkg.add("FiniteDifferences")

# ╔═╡ 17a7cd01-f98e-4cb8-a93b-f80c713d5423
using LinearMaps, LinearAlgebra, MLUtils, LossFunctions

# ╔═╡ 546201a6-096a-4645-9b13-4f563f54a617
using FiniteDifferences

# ╔═╡ 1fc169df-c27d-49c9-9ef8-898289b13ab9
using SparseArrays

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
	if CUDA.functional()
			GeoPhyInv.@init_parallel_stencil(2, true, Float32, 2)
else
	GeoPhyInv.@init_parallel_stencil(2, false, Float32, 2)

end
   
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
medium=Medium(:pizza, 40)

# ╔═╡ aafae9ef-f42c-4a93-a2b5-3de6c59ce8bf
plot(medium)

# ╔═╡ d03a303f-33ef-4d85-8162-092478d512f3
medium0=Medium(:acou_homo2D, 40)

# ╔═╡ 0f2d0963-b91a-4b1c-92d5-31f077e28516
ageom = AGeom(medium.mgrid, :surf, SSrcs(1), Recs(100));

# ╔═╡ d72f7155-b766-46ad-a3af-8390b3cfa1f0
plot(medium0, ageom)

# ╔═╡ 120e51bc-cfd4-4a2a-bcd1-fe942ef224d7
begin
	wav, tgrid = ricker(medium, 2, 2)
	srcwav = SrcWav(tgrid, ageom, [:p]); 
	update!(srcwav, [:p], wav)
end

# ╔═╡ a18fbf29-803c-4f0d-9060-fa1201990b16
plot(srcwav[1])

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

# ╔═╡ d4d425ea-a1f7-42d1-8c22-562cb7260f98
begin
	update!(pa_acoustic, medium0, medium);
	update!(pa_acoustic, pa_acoustic.c.srcwav, [2, 2], verbose=false)
	# # update!(pa_acoustic, medium);
	update!(pa_acoustic)
	dobs = copy(vec(pa_acoustic[:data, 1][:p]));
	xtrue = collect((Iterators.flatten(Array.(pa_acoustic.c.δmod))))
	# , activepw=[1,2], src_flags=[true, false], rec_flags=[false, true]);
end

# ╔═╡ 14e0469e-79f7-4fbe-ba81-b1203d307f73
heatmap(reshape(dobs, :, 100), size=(200, 300))

# ╔═╡ ab2903f4-b077-4cca-8f5c-ad13ae68b7b0
heatmap(Array(pa_acoustic.c.δmod[:K]))

# ╔═╡ 75825692-7105-4382-b099-7736eea44af1
heatmap(pa_acoustic.c.gradients[:K])

# ╔═╡ d23a59fe-45cb-4898-9e63-d4a3c61b0992
heatmap(Array(pa_acoustic.c.δmod[:rhoI]))

# ╔═╡ 6e59506f-9cc7-4db4-9c6c-8b4a38499074
F = LinearMap(pa_acoustic)  	

# ╔═╡ 6bfc9422-b901-4b2b-aab8-0d51af6162e6
GeoPhyInv.Utils.linearmap_test(F)

# ╔═╡ c1839cb1-b9de-4738-8dec-22b74d562934
GeoPhyInv.Utils.adj_test

# ╔═╡ 2b30489a-d03b-4c3b-a75b-3a764156032e
function adjtest()
	x= zeros(size(F, 2)); x[div(length(x), 4)]=rand(); 
	x2= zeros(size(F, 2)); x2[div(length(x), 5)]=rand(); y = F * x2;
	yy = F * x;
	xx = F' * y
	
	(; A=dot(yy, y), B=dot(xx, x))
end

# ╔═╡ 8d92f37d-ee99-4022-8a0f-6006184eb1a0
adjtest()

# ╔═╡ ebe5ecca-dcf6-4bbe-b293-ea374941f844
function update_xs!(x, xs)
  xρ, xinvμ = chunk(x, 2)
  N = length(xρ)
  N2 = div(N, 2)
  X = view(xρ, N2:N2+length(xs)-1)
  copyto!(X, xs)
  return xs
end

# ╔═╡ 12e2dd8d-8422-45da-937f-19b913bb8e0d
begin
	xtest = zeros(size(F, 2));
	update_xs!(xtest, randn(10).*100);
end

# ╔═╡ 9af5dbd9-07c2-4e6c-bff9-3a5dbdf90382
x1 = F' * F * xtest

# ╔═╡ d7e5cee7-20f0-4ffc-b0ec-2bd28317ef29
function get_xs(x, xs)
  xρ, xinvμ = chunk(x, 2)
  N = length(xρ)
  N2 = div(N, 2)
  return xρ[N2:N2+length(xs)-1]
end

# ╔═╡ b11b0764-134b-48af-b90c-5b8ed1f5bc5c
loss = L2DistLoss()

# ╔═╡ 4ecb67e1-3ee6-4da1-bc5d-6c882f704fc7
function Jdata(dobs, d, loss=loss)
  return mapreduce(+, dobs, d) do d1, d2
    value(loss, d1, d2)
  end
end

# ╔═╡ c7187dd0-dfdc-44bc-8657-68db6e33599b
function Js(xs)
	fill!(xtest, 0.0)
	update_xs!(xtest, xs.*1e5)
	yy = F * xtest
	return Jdata(dobs, yy)
end

# ╔═╡ 0006968f-951a-4a05-a191-5e97d2e6f25d
function get_gradient(xs)
		fill!(xtest, 0.0)
	update_xs!(xtest, xs)
d=F*xtest
adj_source = deriv.(loss, dobs, d);
	return get_xs(F' * adj_source, zeros(3))
end

# ╔═╡ 89def4f1-0b75-4fc7-b234-43f5b8194bd8
# g1=grad(central_fdm(2, 1), Js, xs)

# ╔═╡ 63808441-9645-43e5-9f9c-882712d07f6f
xs = get_xs(xtest, zeros(3))

# ╔═╡ f6cc0d3b-591a-4c3e-b202-2fd6cd2767b7
xs

# ╔═╡ e59b245e-fa79-40c3-ba82-b17f1431c562
Js(xs)

# ╔═╡ 23d9fc12-4803-4162-a544-b308744f913f
g2 = get_gradient(xs)

# ╔═╡ ff3e948a-633b-4763-90eb-b51c7a555138
g1[1] ./g2 

# ╔═╡ 0f8685fc-7341-4cd5-a7cb-bf6b5b5d5ab9
nz, nx = length.(pa_acoustic.c.exmedium.mgrid)

# ╔═╡ 158ef323-b376-46d1-aea8-c26061ee5b56
heatmap(reshape(xtrue[1:nx*nz], nz, nx))

# ╔═╡ 3b3ac5f6-4fd0-430e-b35c-84bf0d66c2ea
heatmap(reshape((x1[1:nx*nz]), nz, nx))

# ╔═╡ abeca72d-c414-451a-b119-303d0542d268
reshape(x[1:nx*nz], nz, nx) |> iszero

# ╔═╡ c9a6db62-ee71-4f62-94d0-de43f66e0342


# ╔═╡ 19d90bda-f51e-4c99-8abe-ea5cbe0b39e8
heatmap(reshape(xx[1:nx*nz], nz, nx))

# ╔═╡ f9c12045-7287-4471-a181-50a32625b0cd
func(x) = norm(F*pad(x), 2)^2

# ╔═╡ 3f0b4a46-0742-453a-b8db-0efb9bfb619e
# xs = truncate(x)

# ╔═╡ f9a92623-ade1-4a65-b89b-e9fe9e9cfb22
# gm=grad(forward_fdm(2, 1), func, xs)[1]

# ╔═╡ fca44bc2-7af2-4230-801c-6f574be2b883
plot(gm)# plot!(gm2)

# ╔═╡ ded5d47e-3bde-4355-8628-e1a1b2150c60
gm2=truncate(F'*F*x)

# ╔═╡ 8f364510-de8b-477a-8ab3-68819df1cc68
plot(gm2)

# ╔═╡ 5418320a-9b38-4a52-9b9c-d160fb229c18


# ╔═╡ fd3c3a48-d722-4da4-b558-0599af0d0b1c
heatmap(reshape(y, length(tgrid), :), margin=2Plots.cm)

# ╔═╡ 73ec2fd6-2956-4d92-9a10-d5d7caa456e6
heatmap(reshape(yy, length(tgrid), :), margin=2Plots.cm)

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
# ╠═970c1837-e646-4a3f-9946-9b494af2c21c
# ╠═a34fc72e-8f98-4263-b48e-6498fa69f211
# ╠═78769d02-829e-407a-a211-e0d02928e9a8
# ╠═aafae9ef-f42c-4a93-a2b5-3de6c59ce8bf
# ╠═d03a303f-33ef-4d85-8162-092478d512f3
# ╠═d72f7155-b766-46ad-a3af-8390b3cfa1f0
# ╠═120e51bc-cfd4-4a2a-bcd1-fe942ef224d7
# ╠═0f2d0963-b91a-4b1c-92d5-31f077e28516
# ╠═c54ae16f-7c39-4899-99e6-6d82c8929bde
# ╠═a18fbf29-803c-4f0d-9060-fa1201990b16
# ╠═cc8fddd6-fe22-45d7-8abb-e9c66f7c449c
# ╠═14e0469e-79f7-4fbe-ba81-b1203d307f73
# ╠═d4d425ea-a1f7-42d1-8c22-562cb7260f98
# ╠═ab2903f4-b077-4cca-8f5c-ad13ae68b7b0
# ╠═75825692-7105-4382-b099-7736eea44af1
# ╠═d23a59fe-45cb-4898-9e63-d4a3c61b0992
# ╠═17a7cd01-f98e-4cb8-a93b-f80c713d5423
# ╠═6e59506f-9cc7-4db4-9c6c-8b4a38499074
# ╠═6bfc9422-b901-4b2b-aab8-0d51af6162e6
# ╠═c1839cb1-b9de-4738-8dec-22b74d562934
# ╠═12e2dd8d-8422-45da-937f-19b913bb8e0d
# ╠═158ef323-b376-46d1-aea8-c26061ee5b56
# ╠═9af5dbd9-07c2-4e6c-bff9-3a5dbdf90382
# ╠═f6cc0d3b-591a-4c3e-b202-2fd6cd2767b7
# ╠═3b3ac5f6-4fd0-430e-b35c-84bf0d66c2ea
# ╠═2b30489a-d03b-4c3b-a75b-3a764156032e
# ╠═8d92f37d-ee99-4022-8a0f-6006184eb1a0
# ╠═ebe5ecca-dcf6-4bbe-b293-ea374941f844
# ╠═d7e5cee7-20f0-4ffc-b0ec-2bd28317ef29
# ╠═c7187dd0-dfdc-44bc-8657-68db6e33599b
# ╠═b11b0764-134b-48af-b90c-5b8ed1f5bc5c
# ╠═4ecb67e1-3ee6-4da1-bc5d-6c882f704fc7
# ╠═0006968f-951a-4a05-a191-5e97d2e6f25d
# ╠═e59b245e-fa79-40c3-ba82-b17f1431c562
# ╠═89def4f1-0b75-4fc7-b234-43f5b8194bd8
# ╠═23d9fc12-4803-4162-a544-b308744f913f
# ╠═ff3e948a-633b-4763-90eb-b51c7a555138
# ╠═63808441-9645-43e5-9f9c-882712d07f6f
# ╠═0f8685fc-7341-4cd5-a7cb-bf6b5b5d5ab9
# ╠═abeca72d-c414-451a-b119-303d0542d268
# ╠═c9a6db62-ee71-4f62-94d0-de43f66e0342
# ╠═19d90bda-f51e-4c99-8abe-ea5cbe0b39e8
# ╠═e2e1d696-ad9d-4e87-9e9a-11a3f68d05de
# ╠═546201a6-096a-4645-9b13-4f563f54a617
# ╠═f9c12045-7287-4471-a181-50a32625b0cd
# ╠═1fc169df-c27d-49c9-9ef8-898289b13ab9
# ╠═3f0b4a46-0742-453a-b8db-0efb9bfb619e
# ╠═f9a92623-ade1-4a65-b89b-e9fe9e9cfb22
# ╠═fca44bc2-7af2-4230-801c-6f574be2b883
# ╠═8f364510-de8b-477a-8ab3-68819df1cc68
# ╠═ded5d47e-3bde-4355-8628-e1a1b2150c60
# ╠═5418320a-9b38-4a52-9b9c-d160fb229c18
# ╠═fd3c3a48-d722-4da4-b558-0599af0d0b1c
# ╠═73ec2fd6-2956-4d92-9a10-d5d7caa456e6
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
