### A Pluto.jl notebook ###
# v0.19.21

#> [frontmatter]
#> title = "SeisForwExpt"
#> description = "This notebook demonstrates how an instance of SeisForwExpt can be used iteratively to perform forward modeling"

using Markdown
using InteractiveUtils

# ╔═╡ 5fc7d3bf-418d-4487-a2d6-bc07c736374b
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

# ╔═╡ 98de96ba-8a96-4039-bb9f-f52ad9ee31bc
begin
    Pkg.activate(Base.current_project())
    Pkg.instantiate()
end

# ╔═╡ 9fd95dde-ff9e-490e-b357-d8e0d83823dc
@revise using GeoPhyInv

# ╔═╡ 68d5c87e-02a6-41ed-aac1-559d0bfde12a
TableOfContents()

# ╔═╡ c68a854d-7ae7-4f29-80e3-55207bd61d29
Pkg.activate

# ╔═╡ a6b59ba1-ea71-467f-bf26-49a634a1bbd9
md"""
# Reusing SeisForwExpt
"""

# ╔═╡ 50494fe6-6bdf-4db1-8e61-6c9b39627b1f
md"""
This tutorial demonstrates how an instance of `SeisForwExpt` can be used iteratively 
to perform forward modeling using `update!` for 
various bundles of medium parameters. We will first load a predefined medium, a simple acquisition setting and source wavelets.
"""

# ╔═╡ ae4201a4-7bb1-456c-bdec-d55127324118
# a medium with 5 m spatial sampling
medium = ElasticMedium(Homogeneous(), 5)

# ╔═╡ ded21d43-99d6-453b-a858-e32f3f82e605
# surface seismic acquisition with 3 supersources and 100 receivers
ageom = AGeom(medium.grid, :surf, SSrcs(1), Recs(100))

# ╔═╡ a47b1e5d-a08c-4750-b1f3-f8a62225a19f
# lets choose a time grid
tgrid = range(0.0, stop=2.0, length=2500)

# ╔═╡ dbe71ae1-8053-4b62-b0a8-4aa6810a6655
# load a source wavelet that peaks at 0.25 s and 15 Hz dominant frequency
wav = ricker(15.0, tgrid, tpeak=0.25);

# ╔═╡ 8e94a56e-d483-419c-adb5-7b1889677929
begin
    # allocate an instance of `Srcs` on `:vz` field
    srcwav = Srcs(tgrid, ageom, [:vz])
    # use the Ricker wavelet for all supersources
    update!(srcwav, [:vz], wav)
end

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
    medium=medium,
    ageom=ageom,
    srcwav=srcwav,
    tgrid=tgrid,
    rfields=[:vz],
    verbose=true,
);

# ╔═╡ 9a3360c0-72c7-420f-a7c0-5565a1383f07
md"""
Before modeling, we will now generate a perturbed medium, which is "similar" to the original medium used to create the `Expt` type.
"""

# ╔═╡ df4a4072-a91a-471a-88d9-b74ed58e7823
begin
    medium_box = deepcopy(medium)
    update!(medium_box, [:vp, :rho, :vs], rectangle=[[-500, -500], [500, 500]], perc=5.0) # perturbed velocity box

end

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
    t1 = update!(pa) # run forward simulation
    data = deepcopy(pa[:data]) # copy data 

    update!(pa, medium_box) # load another medium into pa
    t2 = update!(pa) # new simulation with new medium
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
plot(plot(data, 99.9), plot(data_box, 99.9), size=(800, 500))

# ╔═╡ ec7985c6-7072-4ac4-9caa-74a815193750
@test data ≠ data_box

# ╔═╡ Cell order:
# ╠═68d5c87e-02a6-41ed-aac1-559d0bfde12a
# ╠═5fc7d3bf-418d-4487-a2d6-bc07c736374b
# ╠═c68a854d-7ae7-4f29-80e3-55207bd61d29
# ╠═98de96ba-8a96-4039-bb9f-f52ad9ee31bc
# ╠═9fd95dde-ff9e-490e-b357-d8e0d83823dc
# ╟─a6b59ba1-ea71-467f-bf26-49a634a1bbd9
# ╟─50494fe6-6bdf-4db1-8e61-6c9b39627b1f
# ╠═ae4201a4-7bb1-456c-bdec-d55127324118
# ╠═ded21d43-99d6-453b-a858-e32f3f82e605
# ╠═a47b1e5d-a08c-4750-b1f3-f8a62225a19f
# ╠═dbe71ae1-8053-4b62-b0a8-4aa6810a6655
# ╠═8e94a56e-d483-419c-adb5-7b1889677929
# ╠═10be735e-4012-4fa8-b538-9f9a08694458
# ╠═7fa1ddc5-fe10-44e2-b764-533be1840ab6
# ╟─cec954f6-c278-4003-99e4-6c995f74b606
# ╠═c7d96dda-17f5-4dbb-98bf-4c29fd784574
# ╠═9a3360c0-72c7-420f-a7c0-5565a1383f07
# ╠═df4a4072-a91a-471a-88d9-b74ed58e7823
# ╠═07da9113-a518-4925-8978-8407dc60c605
# ╟─7ac4f51d-73b4-4f9c-88a4-a48911742282
# ╠═033cd313-4e87-4364-a1f8-136e0e06bb65
# ╠═fd11d02d-80cd-4373-a8d1-e0d50d93dcf8
# ╠═368c2bbd-94fb-4826-96eb-8fa37a7e35a3
# ╠═7236dc4a-fe04-492d-9532-170ac7c22fb8
# ╠═4612405a-530b-4eab-8b84-40dc907929e4
# ╠═ec7985c6-7072-4ac4-9caa-74a815193750
