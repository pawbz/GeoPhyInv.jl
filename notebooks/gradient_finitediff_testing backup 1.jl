### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ 75ad24f7-8b57-401b-8170-d3e35f3c960f
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

# ╔═╡ 893cd92d-3d0d-4fb2-b465-f61a237fc154
begin
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
end

# ╔═╡ f86410b5-87f7-4761-91a8-856406d00234
@revise using GeoPhyInv

# ╔═╡ bd68d9b1-633b-4c9b-b759-791581714306
using Statistics, LossFunctions, LinearAlgebra, MLUtils, FiniteDifferences, SparseArrays

# ╔═╡ bbfa8c8e-1ef5-459b-b5f5-edb15ff38fa8
TableOfContents()

# ╔═╡ e7842e79-a155-4d0b-b594-1503b2cca8a1
GeoPhyInv._fd_use_gpu

# ╔═╡ de9a4882-cba9-4707-a900-3d82b9a2faa1
# testing parameterization
begin
	pa_mod = GeoPhyInv.SeisForwExpt(:acou_homo2D, npw=2);
	m1=GeoPhyInv.get_modelvector(pa_mod); 
	update!(pa_mod, m1)
	m2=GeoPhyInv.get_modelvector(pa_mod)
	@test m1 ≈ m2
end

# ╔═╡ 15654987-9979-46cb-bae6-1fd631bae58b
Records

# ╔═╡ c6b19c3d-6ccb-436f-a191-d3cc1958a91a
begin
	pa_true = GeoPhyInv.SeisForwExpt(:acou_homo2D);
	update!(pa_true); dobs=deepcopy(pa_true.c.data[1]);
end

# ╔═╡ 43fb3211-9ffc-4ea5-a459-539cf55ba009
pa_inv = SeisInvExpt(pa_mod, dobs, [range(-500, 500, length=5), range(-500, 500, length=5)])

# ╔═╡ 16727899-e02c-4851-bf72-5ceab0078d6c
m=GeoPhyInv.get_modelvector(pa_inv)

# ╔═╡ 4ce63c81-9367-4635-89af-27219f0fe478
loss=L2DistLoss()

# ╔═╡ 57bc5dd6-a33e-4439-ad27-fdbb9b3b17e3
GeoPhyInv.lossvalue(loss, pa_true.c.data[1], pa_mod.c.data[1])

# ╔═╡ 6c037136-4822-438f-85ba-9d61c724f0d9
GeoPhyInv.lossvalue(m, loss, dobs, pa_mod)

# ╔═╡ fc739bbe-7cf6-4701-9462-7813768a924a
plot(dobs[1])

# ╔═╡ 74401c84-0f96-46ff-bbfd-7414440c907b
extrema(dobs[1].d[1])

# ╔═╡ 1dcee841-0a03-422c-874e-ecba38ab4cab
plot(pa_mod.c.srcwav[2][1])

# ╔═╡ c3cc4fd8-7477-47a8-b1ce-3735a1868666
plot(pa_mod.c.data[1][1])

# ╔═╡ f059068e-e5ad-44b7-94b7-c99e1fafaff9
dobs[1].d[:vz] |> extrema

# ╔═╡ f29904f2-c1bc-48e7-9dd8-c58a8c550f16
@show m

# ╔═╡ 8a47e8ba-08df-40bc-ad77-539c8da17943
begin
	g=similar(m);  GeoPhyInv.gradient!(g, m, pa_inv);		
end

# ╔═╡ 4bbcd3e5-62f5-408d-a4a4-eee04fb61b26
gKI = reshape(first(chunk(g, 2)), 5, 5)

# ╔═╡ aa0f7f1b-65ff-4a0b-9fb0-e6660911fdc8
heatmap(Array(deepcopy(pa_mod.c.gradients[:KI])));

# ╔═╡ 630faaff-f7f0-46d7-a36f-31a54e589c8e
heatmap(Array(gKI), c=:seismic)

# ╔═╡ 8ec00fb7-fe18-4741-b987-813ec97d73fe
heatmap(Array(reshape(first(chunk(gm, 2)), 5, 5)), c=:seismic)

# ╔═╡ be01b30c-3918-4c9f-8e9c-f971f9875ace
pa_inv.paf.c.attrib_mod.mode

# ╔═╡ e42f473e-1d05-407a-a7df-3d406ab715c1
heatmap(Array(pa_mod.c.gradients[1]))

# ╔═╡ 417388ad-8709-4d4a-9d81-d4bbb0351a92
heatmap(Array(pa_mod.c.gradients[1]), clim=(-1e-11, 1e-11), c=:seismic)

# ╔═╡ 4bb2c6c1-01ad-44df-9382-635b05b589d0
extrema(pa_mod.c.gradients[1])

# ╔═╡ 4b49262a-28a5-4182-bc46-43e567056bd6
xs=zeros(50)

# ╔═╡ de9d8725-cf72-4913-bff8-b90febbfecdf
xs

# ╔═╡ adee44a1-fb22-4714-9fa4-32f373eb5bf7
eps(Float32)

# ╔═╡ 8b9bb7bf-2c36-41ad-bf8e-5a251e24ec96
function Js(m)
	@show mean(abs.(m))
	return GeoPhyInv.lossvalue(GeoPhyInv.Data.Array(m), pa_inv)
end

# ╔═╡ 22b5fc09-5697-4599-9d27-a64a3584cb0c
Js(GeoPhyInv.Data.Array(xs))

# ╔═╡ c6626e0f-4bf1-4e29-a0c3-1e7f23805f61
# gm=grad(central_fdm(2, 1, factor=1e12), Js, xs)[1]

# ╔═╡ 9371b452-4ad6-4d07-ae98-c86239b6c110
gs=get_xs(g, xs)

# ╔═╡ 19b4cdf4-c87d-42f7-b853-5468b6e2ef52
Array(gm) ./ Array(g)

# ╔═╡ 3a08fc58-311c-4b42-8150-865686b3aeed
step(pa_mod.c.srcwav[1][1].grid)

# ╔═╡ 849b8899-674a-49d8-ad7a-1583b7f9261a
step(pa_mod.c.srcwav[1][1].grid) ./ pa_mod.c.medium[:rho]

# ╔═╡ Cell order:
# ╠═bbfa8c8e-1ef5-459b-b5f5-edb15ff38fa8
# ╠═75ad24f7-8b57-401b-8170-d3e35f3c960f
# ╠═893cd92d-3d0d-4fb2-b465-f61a237fc154
# ╠═f86410b5-87f7-4761-91a8-856406d00234
# ╠═bd68d9b1-633b-4c9b-b759-791581714306
# ╠═e7842e79-a155-4d0b-b594-1503b2cca8a1
# ╠═de9a4882-cba9-4707-a900-3d82b9a2faa1
# ╠═15654987-9979-46cb-bae6-1fd631bae58b
# ╠═43fb3211-9ffc-4ea5-a459-539cf55ba009
# ╠═16727899-e02c-4851-bf72-5ceab0078d6c
# ╠═c6b19c3d-6ccb-436f-a191-d3cc1958a91a
# ╠═4ce63c81-9367-4635-89af-27219f0fe478
# ╠═57bc5dd6-a33e-4439-ad27-fdbb9b3b17e3
# ╠═09854dd9-bc7d-4eca-8f33-0d89f77ae5d8
# ╠═6c037136-4822-438f-85ba-9d61c724f0d9
# ╠═fc739bbe-7cf6-4701-9462-7813768a924a
# ╠═74401c84-0f96-46ff-bbfd-7414440c907b
# ╠═1dcee841-0a03-422c-874e-ecba38ab4cab
# ╠═c3cc4fd8-7477-47a8-b1ce-3735a1868666
# ╠═f059068e-e5ad-44b7-94b7-c99e1fafaff9
# ╠═f29904f2-c1bc-48e7-9dd8-c58a8c550f16
# ╠═8a47e8ba-08df-40bc-ad77-539c8da17943
# ╠═4bbcd3e5-62f5-408d-a4a4-eee04fb61b26
# ╠═aa0f7f1b-65ff-4a0b-9fb0-e6660911fdc8
# ╠═630faaff-f7f0-46d7-a36f-31a54e589c8e
# ╠═8ec00fb7-fe18-4741-b987-813ec97d73fe
# ╠═be01b30c-3918-4c9f-8e9c-f971f9875ace
# ╠═e42f473e-1d05-407a-a7df-3d406ab715c1
# ╠═417388ad-8709-4d4a-9d81-d4bbb0351a92
# ╠═4bb2c6c1-01ad-44df-9382-635b05b589d0
# ╠═4b49262a-28a5-4182-bc46-43e567056bd6
# ╠═de9d8725-cf72-4913-bff8-b90febbfecdf
# ╠═adee44a1-fb22-4714-9fa4-32f373eb5bf7
# ╠═22b5fc09-5697-4599-9d27-a64a3584cb0c
# ╠═8b9bb7bf-2c36-41ad-bf8e-5a251e24ec96
# ╠═c6626e0f-4bf1-4e29-a0c3-1e7f23805f61
# ╠═9371b452-4ad6-4d07-ae98-c86239b6c110
# ╠═19b4cdf4-c87d-42f7-b853-5468b6e2ef52
# ╠═3a08fc58-311c-4b42-8150-865686b3aeed
# ╠═849b8899-674a-49d8-ad7a-1583b7f9261a
