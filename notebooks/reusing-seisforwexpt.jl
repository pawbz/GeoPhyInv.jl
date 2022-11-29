### A Pluto.jl notebook ###
# v0.19.15

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
using PlutoLinks: @revise
using PlutoUI, PlutoTest
end

# ╔═╡ 981f55af-1557-49c5-921d-2e7e343a511b
using Plots, Statistics

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

# ╔═╡ de2f97fa-f64d-4754-bd34-e44dbf13c336
md"""
In order to install `GeoPhyInv` enter these package manager commands in the REPL.
```julia
using Pkg
Pkg.add(PackageSpec(name="GeoPhyInv",url="https://github.com/pawbz/GeoPhyInv.jl.git"))
```
It is important to configure GeoPhyInv with a macro `@init_parallel_stencil` before anything else. If you need to change this configuration, the julia kernel must be restarted.
```julia
using GeoPhyInv; @init_parallel_stencil(⋯)
```
"""

# ╔═╡ d9b71485-8b64-4ad0-a242-dde6300af835
# ╠═╡ show_logs = false
begin
	reload_geophyinv
	@init_parallel_stencil(2, false, Float32, 2)
end

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

# ╔═╡ 5e13eef1-e6c9-4a6a-ac47-62ecac5ba67b
Medium

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
begin
	# allocate an instance of `SrcWav` on `:vz` field
	srcwav = SrcWav(tgrid, ageom, [:vz]); 
	# use the Ricker wavelet for all supersources
	update!(srcwav, [:vz], wav);
end

# ╔═╡ d4aab91a-74da-450e-b610-66858a0ef33c


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
plot(plot(data, 99.9), plot(data_box, 99.9), size=(800, 500))

# ╔═╡ ec7985c6-7072-4ac4-9caa-74a815193750
@test data ≠ data_box

# ╔═╡ Cell order:
# ╟─86d3f068-a979-42f5-a9e7-138e94c16b38
# ╟─de2f97fa-f64d-4754-bd34-e44dbf13c336
# ╟─0ffd8ee4-1735-4756-befb-c7c10d08eb34
# ╟─7687d367-f9d5-4539-9156-e26d87379f87
# ╟─d9b71485-8b64-4ad0-a242-dde6300af835
# ╠═981f55af-1557-49c5-921d-2e7e343a511b
# ╟─a6b59ba1-ea71-467f-bf26-49a634a1bbd9
# ╟─50494fe6-6bdf-4db1-8e61-6c9b39627b1f
# ╠═5e13eef1-e6c9-4a6a-ac47-62ecac5ba67b
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
# ╠═ec7985c6-7072-4ac4-9caa-74a815193750
