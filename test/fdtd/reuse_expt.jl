### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 98bb6ef8-dff3-4b53-a2bb-7cfb862fe098
import Pkg; Pkg.activate("GeoPhyInv")

# ╔═╡ 08eb078a-0c2a-11ed-2120-39fc17ba6f77
begin
	using Revise
	using GeoPhyInv
	using Statistics
	using Test
	using CUDA
	using PlutoUI
	CUDA.device!(1)
	@init_parallel_stencil(2, true, Float32, 2)
end

# ╔═╡ 869c583c-c24f-4451-a3c4-056963e88445
using Plots; gr()

# ╔═╡ 50494fe6-6bdf-4db1-8e61-6c9b39627b1f
md"""
This tutorial demonstrates how an instance of `SeisForwExpt` can be used iteratively 
to perform forward modeling using `update!` for 
various bundles of medium parameters. We will first load a predefined medium, a simple acquisition setting and source wavelets.
"""

# ╔═╡ ae4201a4-7bb1-456c-bdec-d55127324118
# a medium with 5 m spatial sampling
medium = Medium(:elastic_homo2D, 5)

# ╔═╡ ded21d43-99d6-453b-a858-e32f3f82e605
# surface seismic acquisition with 3 supersources and 100 receivers
ageom = AGeom(medium.mgrid, :surf, SSrcs(3), Recs(100))

# ╔═╡ a47b1e5d-a08c-4750-b1f3-f8a62225a19f
# lets choose a time grid
tgrid = range(0.0, stop = 2.0, length = 2500)

# ╔═╡ dbe71ae1-8053-4b62-b0a8-4aa6810a6655
# load a source wavelet that peaks at 0.25 s and 15 Hz dominant frequency
wav = ricker(15.0, tgrid, tpeak = 0.25); 

# ╔═╡ 8e94a56e-d483-419c-adb5-7b1889677929
# allocate an instance of `SrcWav` on `:vz` field
srcwav = SrcWav(tgrid, ageom, [:vz]); 

# ╔═╡ d4aab91a-74da-450e-b610-66858a0ef33c
# use the Ricker wavelet for all supersources
update!(srcwav, [:vz], wav);

# ╔═╡ 10be735e-4012-4fa8-b538-9f9a08694458
# lets plot the medium
plot(medium, ageom)

# ╔═╡ 7fa1ddc5-fe10-44e2-b764-533be1840ab6
# plot source wavelet and its spectrum
plot(srcwav[1])

# ╔═╡ cec954f6-c278-4003-99e4-6c995f74b606
md"""
Now we have all the stuff necessary to allocate an instance of `SeisForwExpt`. 
"""

# ╔═╡ c7d96dda-17f5-4dbb-98bf-4c29fd784574
pa = SeisForwExpt(
    FdtdElastic(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    tgrid = tgrid,
    rfields = [:vz],
    verbose = true,
); 

# ╔═╡ 9a3360c0-72c7-420f-a7c0-5565a1383f07
md"""
Before modeling, we will now generate a perturbed medium, which is "similar" to the original medium used to create the `Expt` type.
"""

# ╔═╡ df4a4072-a91a-471a-88d9-b74ed58e7823
medium_box = similar(medium); update!(medium_box, [:vp, :rho, :vs], rectangle = [[-500, -500], [500, 500]], perc = 5.0) # perturbed velocity box


# ╔═╡ 07da9113-a518-4925-8978-8407dc60c605
# lets plot it
plot(medium_box, ageom)

# ╔═╡ 7ac4f51d-73b4-4f9c-88a4-a48911742282
md"""
Finally, we will now run two simulations in different media using the same `Expt` instance. It is possible to change the medium in `pa` without memory allocation.
"""

# ╔═╡ 033cd313-4e87-4364-a1f8-136e0e06bb65
begin
	update!(pa, medium) # load medium into pa
	t1=update!(pa) # run forward simulation
	data = deepcopy(pa[:data]) # copy data 

	update!(pa, medium_box) # load another medium into pa
	t2=update!(pa) # new simulation with new medium
	data_box = deepcopy(pa[:data]) # copy data
end

# ╔═╡ fd11d02d-80cd-4373-a8d1-e0d50d93dcf8
md"""
Lets look at the timer objects.
"""

# ╔═╡ 368c2bbd-94fb-4826-96eb-8fa37a7e35a3
t1

# ╔═╡ 7236dc4a-fe04-492d-9532-170ac7c22fb8
t2

# ╔═╡ 4612405a-530b-4eab-8b84-40dc907929e4
plot(data, 99)

# ╔═╡ f3d41a34-1ef4-45ac-ab74-aed418fb711a
plot(data_box, 99)

# ╔═╡ ec7985c6-7072-4ac4-9caa-74a815193750
@test data ≠ data_box

# ╔═╡ b91b6889-6d13-478d-952a-e1beef46b867
md"""
Notebook configuration.
"""

# ╔═╡ 634c9676-53e9-4fa0-bfa6-2f5250b5964f


# ╔═╡ Cell order:
# ╠═50494fe6-6bdf-4db1-8e61-6c9b39627b1f
# ╠═ae4201a4-7bb1-456c-bdec-d55127324118
# ╠═ded21d43-99d6-453b-a858-e32f3f82e605
# ╠═a47b1e5d-a08c-4750-b1f3-f8a62225a19f
# ╠═dbe71ae1-8053-4b62-b0a8-4aa6810a6655
# ╠═8e94a56e-d483-419c-adb5-7b1889677929
# ╠═d4aab91a-74da-450e-b610-66858a0ef33c
# ╠═10be735e-4012-4fa8-b538-9f9a08694458
# ╠═7fa1ddc5-fe10-44e2-b764-533be1840ab6
# ╠═cec954f6-c278-4003-99e4-6c995f74b606
# ╠═c7d96dda-17f5-4dbb-98bf-4c29fd784574
# ╠═9a3360c0-72c7-420f-a7c0-5565a1383f07
# ╠═df4a4072-a91a-471a-88d9-b74ed58e7823
# ╠═07da9113-a518-4925-8978-8407dc60c605
# ╠═7ac4f51d-73b4-4f9c-88a4-a48911742282
# ╠═033cd313-4e87-4364-a1f8-136e0e06bb65
# ╠═fd11d02d-80cd-4373-a8d1-e0d50d93dcf8
# ╠═368c2bbd-94fb-4826-96eb-8fa37a7e35a3
# ╠═7236dc4a-fe04-492d-9532-170ac7c22fb8
# ╠═4612405a-530b-4eab-8b84-40dc907929e4
# ╠═f3d41a34-1ef4-45ac-ab74-aed418fb711a
# ╠═ec7985c6-7072-4ac4-9caa-74a815193750
# ╠═b91b6889-6d13-478d-952a-e1beef46b867
# ╠═869c583c-c24f-4451-a3c4-056963e88445
# ╠═08eb078a-0c2a-11ed-2120-39fc17ba6f77
# ╠═634c9676-53e9-4fa0-bfa6-2f5250b5964f
# ╠═98bb6ef8-dff3-4b53-a2bb-7cfb862fe098
