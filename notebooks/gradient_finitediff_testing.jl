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

# ╔═╡ de9a4882-cba9-4707-a900-3d82b9a2faa1
# testing parameterization
begin
	pa_mod = GeoPhyInv.SeisForwExpt(:acou_homo2D, npw=2);
	# m1=GeoPhyInv.get_modelvector(pa_mod, [:KI]); 
	# @time update!(pa_mod, m1, [:KI])
	# m2=GeoPhyInv.get_modelvector(pa_mod, [:KI])
	# @test m1 ≈ m2
end

# ╔═╡ 2a737f8b-92d0-4042-96eb-6b85ea4030c7
# @time update!(m1, pa_mod, [:KI])

# ╔═╡ 9d3aeb6f-936d-4173-b41c-2159a458b8ce
m1= GeoPhyInv.get_modelvector(pa_mod, [:KI])

# ╔═╡ 0ecf589f-a5e1-4139-8ec4-563f4a7c38f7


# ╔═╡ f49439df-3db3-48ae-bf9d-1387e3c14c8f
N=2

# ╔═╡ c6b19c3d-6ccb-436f-a191-d3cc1958a91a
begin
	pa_true = GeoPhyInv.SeisForwExpt(:acou_homo2D);
	update!(pa_true); dobs=deepcopy(pa_true.c.data[1]);
end

# ╔═╡ 43fb3211-9ffc-4ea5-a459-539cf55ba009
pa_inv = SeisInvExpt(pa_mod, dobs, [range(-50, 50, length=N), range(-50, 50, length=N)], [:KI])

# ╔═╡ 4ce63c81-9367-4635-89af-27219f0fe478
loss=L2DistLoss()

# ╔═╡ 57bc5dd6-a33e-4439-ad27-fdbb9b3b17e3
GeoPhyInv.lossvalue(loss, pa_true.c.data[1], pa_mod.c.data[1])

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

# ╔═╡ 6e7d2876-4345-441d-9d29-e36fa02ca2eb
GeoPhyInv.get_modelvector(pa_inv)

# ╔═╡ aa0f7f1b-65ff-4a0b-9fb0-e6660911fdc8
heatmap(Array(deepcopy(pa_mod.c.gradients[:KI])));

# ╔═╡ 839cf146-639c-4d13-9470-d25882837169
1/prod(step.(pa_true.c.medium.mgrid)) * pa_true.c.srcwav[1][1].grid |> step

# ╔═╡ d3b4cd08-eca1-4fac-915f-a64dc29b870f
pa_true.c.srcwav[1][1].grid |> step

# ╔═╡ 911518e0-0fc1-43fe-8df9-b4e27382181f


# ╔═╡ 417388ad-8709-4d4a-9d81-d4bbb0351a92
heatmap(Array(pa_mod.c.gradients[1]), clim=(-1e-11, 1e-11), c=:seismic)

# ╔═╡ 4b49262a-28a5-4182-bc46-43e567056bd6
xs=zeros(N*N)

# ╔═╡ c895f6d7-2012-4929-85dd-f0fa1bb97866
m=GeoPhyInv.Data.Array(xs)

# ╔═╡ 6c037136-4822-438f-85ba-9d61c724f0d9
GeoPhyInv.lossvalue(m, loss, dobs, pa_mod, [:KI])

# ╔═╡ b57d5882-660b-4b71-b31e-8d461d70d14f
m

# ╔═╡ 8a47e8ba-08df-40bc-ad77-539c8da17943
begin
	g=similar(m);  GeoPhyInv.gradient!(g, m, pa_inv);		
	heatmap(Array(pa_mod.c.gradients[1]))
end

# ╔═╡ cd5bff59-6466-43d7-bd10-795d6f55d4ef
heatmap(Array(reshape(g, N, N)), c=:seismic)

# ╔═╡ 4bbcd3e5-62f5-408d-a4a4-eee04fb61b26
gKI = reshape(g, N, N)

# ╔═╡ 630faaff-f7f0-46d7-a36f-31a54e589c8e
heatmap(Array(gKI), c=:seismic)

# ╔═╡ 8b9bb7bf-2c36-41ad-bf8e-5a251e24ec96
function Js(m)
	return GeoPhyInv.lossvalue(GeoPhyInv.Data.Array(m), pa_inv)
end

# ╔═╡ 22b5fc09-5697-4599-9d27-a64a3584cb0c
Js(GeoPhyInv.Data.Array(xs))

# ╔═╡ c6626e0f-4bf1-4e29-a0c3-1e7f23805f61
gm=grad(central_fdm(2, 1, factor=1e7), Js, xs)[1]

# ╔═╡ 8ec00fb7-fe18-4741-b987-813ec97d73fe
heatmap(Array(reshape(gm, N, N)), c=:seismic)

# ╔═╡ c6007b1a-3819-4543-bab4-f0dc0c12e075
gm1 = Array(reshape(gm, N, N))

# ╔═╡ 60a9adc0-466d-4e54-9838-f01b0a2c50c5
Array(gKI) ./ gm1

# ╔═╡ 9371b452-4ad6-4d07-ae98-c86239b6c110
gs=get_xs(g, xs)

# ╔═╡ 48356c93-888b-442e-b61a-5f25f24c3f8b
pa_mod.c.ref_mod

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
# ╠═de9a4882-cba9-4707-a900-3d82b9a2faa1
# ╠═2a737f8b-92d0-4042-96eb-6b85ea4030c7
# ╠═9d3aeb6f-936d-4173-b41c-2159a458b8ce
# ╠═0ecf589f-a5e1-4139-8ec4-563f4a7c38f7
# ╠═f49439df-3db3-48ae-bf9d-1387e3c14c8f
# ╠═43fb3211-9ffc-4ea5-a459-539cf55ba009
# ╠═c6b19c3d-6ccb-436f-a191-d3cc1958a91a
# ╠═4ce63c81-9367-4635-89af-27219f0fe478
# ╠═57bc5dd6-a33e-4439-ad27-fdbb9b3b17e3
# ╠═6c037136-4822-438f-85ba-9d61c724f0d9
# ╠═fc739bbe-7cf6-4701-9462-7813768a924a
# ╠═74401c84-0f96-46ff-bbfd-7414440c907b
# ╠═1dcee841-0a03-422c-874e-ecba38ab4cab
# ╠═c3cc4fd8-7477-47a8-b1ce-3735a1868666
# ╠═f059068e-e5ad-44b7-94b7-c99e1fafaff9
# ╠═6e7d2876-4345-441d-9d29-e36fa02ca2eb
# ╠═c895f6d7-2012-4929-85dd-f0fa1bb97866
# ╠═b57d5882-660b-4b71-b31e-8d461d70d14f
# ╠═8a47e8ba-08df-40bc-ad77-539c8da17943
# ╠═aa0f7f1b-65ff-4a0b-9fb0-e6660911fdc8
# ╠═630faaff-f7f0-46d7-a36f-31a54e589c8e
# ╠═8ec00fb7-fe18-4741-b987-813ec97d73fe
# ╠═cd5bff59-6466-43d7-bd10-795d6f55d4ef
# ╠═4bbcd3e5-62f5-408d-a4a4-eee04fb61b26
# ╠═c6007b1a-3819-4543-bab4-f0dc0c12e075
# ╠═839cf146-639c-4d13-9470-d25882837169
# ╠═d3b4cd08-eca1-4fac-915f-a64dc29b870f
# ╠═911518e0-0fc1-43fe-8df9-b4e27382181f
# ╠═60a9adc0-466d-4e54-9838-f01b0a2c50c5
# ╠═417388ad-8709-4d4a-9d81-d4bbb0351a92
# ╠═4b49262a-28a5-4182-bc46-43e567056bd6
# ╠═22b5fc09-5697-4599-9d27-a64a3584cb0c
# ╠═8b9bb7bf-2c36-41ad-bf8e-5a251e24ec96
# ╠═c6626e0f-4bf1-4e29-a0c3-1e7f23805f61
# ╠═9371b452-4ad6-4d07-ae98-c86239b6c110
# ╠═48356c93-888b-442e-b61a-5f25f24c3f8b
# ╠═19b4cdf4-c87d-42f7-b853-5468b6e2ef52
# ╠═3a08fc58-311c-4b42-8150-865686b3aeed
# ╠═849b8899-674a-49d8-ad7a-1583b7f9261a
