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
    pa_mod = SeisForwExpt(FdtdAcoustic{FullWave}(:forward, 2), Homogeneous())
    # m1=GeoPhyInv.get_modelvector(pa_mod, [:invK]); 
    # @time update!(pa_mod, m1, [:invK])
    # m2=GeoPhyInv.get_modelvector(pa_mod, [:invK])
    # @test m1 ≈ m2
end

# ╔═╡ 2a737f8b-92d0-4042-96eb-6b85ea4030c7
# @time update!(m1, pa_mod, [:invK])

# ╔═╡ f49439df-3db3-48ae-bf9d-1387e3c14c8f
N = 200

# ╔═╡ 6d47ae1d-2554-4d42-af6c-02ce952bf14e
@bind mpara Select([:invK, :rho])

# ╔═╡ 9d3aeb6f-936d-4173-b41c-2159a458b8ce
m1 = GeoPhyInv.get_modelvector(pa_mod, [mpara])

# ╔═╡ 5582110f-327c-419a-8759-ca11844c0c96
mpara

# ╔═╡ 43fb3211-9ffc-4ea5-a459-539cf55ba009
# pa_inv = SeisInvExpt(pa_mod, dobs, [N, N], [mpara])

# ╔═╡ 5297741a-60dc-45b8-8334-bbbc9a2a96c7
aa=randn(5,4)

# ╔═╡ 802fb101-09f6-44b6-9a46-4fa94724b6f1
vv=[1, 0, 0, 0]

# ╔═╡ 110ee9ed-6e28-4b15-8102-a75e42d26d07
stack([randn(3), randn(3)], dims=2)

# ╔═╡ 8e9a92dc-dc58-40d7-9a38-df08852a33a0
typeof(pa_mod.c.medium.grid)

# ╔═╡ 9706aa59-0286-43cf-bdb7-a3d040e0e2bc


# ╔═╡ 247896d7-6c39-4bf6-856b-7afc6f4d6655
supertype(StepRangeLen)

# ╔═╡ c6b19c3d-6ccb-436f-a191-d3cc1958a91a
begin
    pa_true = SeisForwExpt(FdtdAcoustic{FullWave}(:forward, 1), RandScatterer())
    update!(pa_true)
    dobs = deepcopy(pa_true.c.data[1])
end

# ╔═╡ 26c9f7a1-f18b-4fe3-a607-ad1f6a904bb8
pa_inv, pa_src = SeisInvExpt(pa_mod, dobs, [range(-100,100,length=N), range(-100,100,length=N)], [mpara])

# ╔═╡ 0e5eb59a-aad6-4c84-a0e4-c8558964719f
update!(pa_inv, pa_src)

# ╔═╡ 88481dbd-7e07-4b7b-b849-8afa360ed0c3
plot(pa_src[1].s)

# ╔═╡ 4d0ddfc2-02a4-4ded-ab9e-61be16204da7
heatmap(pa_src[1].d)

# ╔═╡ 0c2679ff-a90c-4fec-b349-3872ce8110a8
heatmap(pa_src[1].g)

# ╔═╡ fa46977d-9574-40cb-a433-4e074c39eb1e
heatmap(pa_inv.dobs[1].d[1] - pa_inv.dobs0[1].d[1])

# ╔═╡ 383aa4a2-df34-49c4-9118-2d3ab367a31b
heatmap(pa_inv.paf.c.data[1][1].d[1])

# ╔═╡ 530c8c00-f067-47ba-891a-6ba618c69d22
pa_inv.paf.c.data[1][1].d.array

# ╔═╡ e6339233-2d9f-45d0-bcb1-2320a4ad28fc
dobs

# ╔═╡ b99b6e7a-d037-4c6f-8043-fdb49102503b
pa_true.c.srcwav[1][1].grid

# ╔═╡ 4ce63c81-9367-4635-89af-27219f0fe478
loss = L2DistLoss()

# ╔═╡ 57bc5dd6-a33e-4439-ad27-fdbb9b3b17e3
GeoPhyInv.lossvalue(loss, pa_true.c.data[1], pa_mod.c.data[1])

# ╔═╡ fc739bbe-7cf6-4701-9462-7813768a924a
plot(dobs[1])

# ╔═╡ 8d16ea50-51d9-4d76-9552-ad1a20bfa67b
heatmap(dobs[1][:vz] .- pa_mod.c.data[1][1][:vz])

# ╔═╡ 9530c17a-f19a-4fb1-a045-033f70f54fa4
heatmap(dobs[1][:vz])

# ╔═╡ 74401c84-0f96-46ff-bbfd-7414440c907b
extrema(dobs[1].d[1])

# ╔═╡ 1dcee841-0a03-422c-874e-ecba38ab4cab
plot(pa_mod.c.srcwav[2][1])

# ╔═╡ c3cc4fd8-7477-47a8-b1ce-3735a1868666
plot(pa_mod.c.data[1][1])

# ╔═╡ f059068e-e5ad-44b7-94b7-c99e1fafaff9
dobs[1].d[:vz] |> extrema

# ╔═╡ 6e7d2876-4345-441d-9d29-e36fa02ca2eb
m = GeoPhyInv.get_modelvector(pa_inv)

# ╔═╡ 6c037136-4822-438f-85ba-9d61c724f0d9
GeoPhyInv.lossvalue(m, loss, dobs, pa_mod, [mpara])

# ╔═╡ b57d5882-660b-4b71-b31e-8d461d70d14f
m

# ╔═╡ 8a47e8ba-08df-40bc-ad77-539c8da17943
begin
    g = similar(m)
    GeoPhyInv.gradient!(g, m, pa_inv)
    heatmap(Array(pa_mod.c.gradients[1]))
end

# ╔═╡ aa0f7f1b-65ff-4a0b-9fb0-e6660911fdc8
heatmap(Array(deepcopy(pa_mod.c.gradients[mpara])));

# ╔═╡ 02e4d288-bde3-4ef4-877f-5947ef50816a
heatmap(Array(pa_mod.c.mod[1]))

# ╔═╡ cd4bb6a4-91e4-4ab5-b462-1d70e781032a
Nz, Nx = length.(pa_mod.c.exmedium.grid)

# ╔═╡ cb409839-d520-4d47-8701-48a4ece49990
heatmap(Array(reshape(pa_inv.P' * m, Nz, Nx)))

# ╔═╡ 01552458-5de7-4627-b417-21745db2ade1
heatmap(Array(reshape(m, 200, 200)))

# ╔═╡ 8ec00fb7-fe18-4741-b987-813ec97d73fe
heatmap(Array(reshape(gm, N, N)), c=:seismic)

# ╔═╡ cd5bff59-6466-43d7-bd10-795d6f55d4ef
heatmap(Array(reshape(g, N, N)), c=:seismic, clim=(-1e-15, 1e-15))

# ╔═╡ cb6ba7e9-0685-4e59-96ef-a8a7af640107
mul!

# ╔═╡ 4bbcd3e5-62f5-408d-a4a4-eee04fb61b26
gKI = reshape(g, N, N)

# ╔═╡ 630faaff-f7f0-46d7-a36f-31a54e589c8e
heatmap(Array(gKI), c=:seismic)

# ╔═╡ c6007b1a-3819-4543-bab4-f0dc0c12e075
gm1 = Array(reshape(gm, N, N))

# ╔═╡ 839cf146-639c-4d13-9470-d25882837169
# prod(step.(pa_true.c.medium.grid)) #
pa_true.c.srcwav[1][1].grid |> step |> abs2

# ╔═╡ d3b4cd08-eca1-4fac-915f-a64dc29b870f
GeoPhyInv.Fields(pa_mod.c.attrib_mod, "v")

# ╔═╡ b14a5a36-94c4-4a6c-8b26-730919e3e5ae
filter(x -> x ∉ Fields(pa_mod.c.attrib_mod, "d"), Fields(pa_mod.c.attrib_mod, "v"))

# ╔═╡ 60a9adc0-466d-4e54-9838-f01b0a2c50c5
Array(gKI) ./ gm1

# ╔═╡ 417388ad-8709-4d4a-9d81-d4bbb0351a92
heatmap(Array(pa_mod.c.gradients[1]), clim=(-1e-11, 1e-11), c=:seismic)

# ╔═╡ 4b49262a-28a5-4182-bc46-43e567056bd6
xs = fill(1e-1, N * N)

# ╔═╡ 8b9bb7bf-2c36-41ad-bf8e-5a251e24ec96
function Js(m)
    return GeoPhyInv.lossvalue(GeoPhyInv.Data.Array(m), pa_inv)
end

# ╔═╡ 22b5fc09-5697-4599-9d27-a64a3584cb0c
Js(GeoPhyInv.Data.Array(xs))

# ╔═╡ c6626e0f-4bf1-4e29-a0c3-1e7f23805f61
# gm = grad(central_fdm(2, 1, factor=1e8), Js, xs)[1]

# ╔═╡ 9371b452-4ad6-4d07-ae98-c86239b6c110
gs = get_xs(g, xs)

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
# ╠═f49439df-3db3-48ae-bf9d-1387e3c14c8f
# ╠═6d47ae1d-2554-4d42-af6c-02ce952bf14e
# ╠═5582110f-327c-419a-8759-ca11844c0c96
# ╠═26c9f7a1-f18b-4fe3-a607-ad1f6a904bb8
# ╠═0e5eb59a-aad6-4c84-a0e4-c8558964719f
# ╠═88481dbd-7e07-4b7b-b849-8afa360ed0c3
# ╠═4d0ddfc2-02a4-4ded-ab9e-61be16204da7
# ╠═0c2679ff-a90c-4fec-b349-3872ce8110a8
# ╠═fa46977d-9574-40cb-a433-4e074c39eb1e
# ╠═383aa4a2-df34-49c4-9118-2d3ab367a31b
# ╠═530c8c00-f067-47ba-891a-6ba618c69d22
# ╠═43fb3211-9ffc-4ea5-a459-539cf55ba009
# ╠═5297741a-60dc-45b8-8334-bbbc9a2a96c7
# ╠═802fb101-09f6-44b6-9a46-4fa94724b6f1
# ╠═110ee9ed-6e28-4b15-8102-a75e42d26d07
# ╠═e6339233-2d9f-45d0-bcb1-2320a4ad28fc
# ╠═8e9a92dc-dc58-40d7-9a38-df08852a33a0
# ╠═9706aa59-0286-43cf-bdb7-a3d040e0e2bc
# ╠═247896d7-6c39-4bf6-856b-7afc6f4d6655
# ╠═c6b19c3d-6ccb-436f-a191-d3cc1958a91a
# ╠═b99b6e7a-d037-4c6f-8043-fdb49102503b
# ╠═4ce63c81-9367-4635-89af-27219f0fe478
# ╠═57bc5dd6-a33e-4439-ad27-fdbb9b3b17e3
# ╠═6c037136-4822-438f-85ba-9d61c724f0d9
# ╠═fc739bbe-7cf6-4701-9462-7813768a924a
# ╠═8d16ea50-51d9-4d76-9552-ad1a20bfa67b
# ╠═9530c17a-f19a-4fb1-a045-033f70f54fa4
# ╠═74401c84-0f96-46ff-bbfd-7414440c907b
# ╠═1dcee841-0a03-422c-874e-ecba38ab4cab
# ╠═c3cc4fd8-7477-47a8-b1ce-3735a1868666
# ╠═f059068e-e5ad-44b7-94b7-c99e1fafaff9
# ╠═6e7d2876-4345-441d-9d29-e36fa02ca2eb
# ╠═c895f6d7-2012-4929-85dd-f0fa1bb97866
# ╠═b57d5882-660b-4b71-b31e-8d461d70d14f
# ╠═8a47e8ba-08df-40bc-ad77-539c8da17943
# ╠═aa0f7f1b-65ff-4a0b-9fb0-e6660911fdc8
# ╠═02e4d288-bde3-4ef4-877f-5947ef50816a
# ╠═cd4bb6a4-91e4-4ab5-b462-1d70e781032a
# ╠═cb409839-d520-4d47-8701-48a4ece49990
# ╠═01552458-5de7-4627-b417-21745db2ade1
# ╠═630faaff-f7f0-46d7-a36f-31a54e589c8e
# ╠═8ec00fb7-fe18-4741-b987-813ec97d73fe
# ╠═cd5bff59-6466-43d7-bd10-795d6f55d4ef
# ╠═cb6ba7e9-0685-4e59-96ef-a8a7af640107
# ╠═4bbcd3e5-62f5-408d-a4a4-eee04fb61b26
# ╠═c6007b1a-3819-4543-bab4-f0dc0c12e075
# ╠═839cf146-639c-4d13-9470-d25882837169
# ╠═d3b4cd08-eca1-4fac-915f-a64dc29b870f
# ╠═b14a5a36-94c4-4a6c-8b26-730919e3e5ae
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
