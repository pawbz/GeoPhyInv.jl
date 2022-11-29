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
Pkg.add("Plots")
using PlutoLinks: @revise
using PlutoUI, PlutoTest, Plots
end

# ╔═╡ 981f55af-1557-49c5-921d-2e7e343a511b
using Statistics

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
## Medium
"""

# ╔═╡ 3c549ecd-52dc-47d2-b846-620148089296
@doc Medium

# ╔═╡ addd0bef-5d52-4e06-9350-f6dade3db428
md"## Examples"

# ╔═╡ 50494fe6-6bdf-4db1-8e61-6c9b39627b1f
md"""
Let's load a predefined medium.
"""

# ╔═╡ 80109fc0-7446-47eb-b912-1f767d824fbd
medium = Medium(:elastic_homo3D)

# ╔═╡ 701d36a4-5f09-43ab-bd09-aefc92196182
names(medium)

# ╔═╡ 63ec2c93-539e-4a58-a10a-fb61602ae03e
medium.ref |> println

# ╔═╡ d356c66b-b07e-415a-a763-9d7f235f1c48
medium.bounds |> println

# ╔═╡ bdac654f-5743-477c-af43-17cb8ffc9b60
medium2D = Medium(:elastic_homo2D)

# ╔═╡ 3cc49ee4-06f6-4459-9f42-8838963f344a
plot(medium2D, fields=[:vp, :vs])

# ╔═╡ 15b42a58-5ad7-42a6-9dd4-7db648c33f55
mgrid = [range(0.0, stop = 10.0, step = 0.1), range(0.0, stop = 30.0, step = 0.2)];

# ╔═╡ fc58c853-2665-45f8-997d-1b3217fa8ec0
md"To construct an instance of 2-D Medium, we need ranges for `z` and `x`."

# ╔═╡ 7a721ee1-227b-4975-b4c7-1ac1ec83e973
medium_custom = Medium(mgrid, [:vp, :rho, :vs]);

# ╔═╡ e55ed529-8670-4326-ac26-e3eaa55f624b
begin
	vpb = [2100.0, 2200.0];
	vsb = [1500, 1700];
	rhob = [2100.0, 2300.0];
	update!(medium_custom, [:vp, :vs, :rho], [vpb, vsb, rhob]);
end

# ╔═╡ 566a79e2-2e79-40db-8a46-e0a341989bb0
fill!(medium_custom);

# ╔═╡ bc26144a-167b-4e99-a631-d7e16381fe97
md"Once the basic medium parameters are input, we can access some other derived parameters."

# ╔═╡ 3b9d77ed-ad38-429d-b6de-34722a3850a8
medium[:Zp];

# ╔═╡ 5bf36510-ddb5-43bb-a3a0-02a22578b27a
begin
	medium[:vp] .= 3000.0;
	medium[:vs] .= 2000.0;
	
	println(medium)
end

# ╔═╡ 315f5412-74ec-4b18-9286-9201c7ab5dc1
update!(medium2D, [:vp, :rho], randn_perc = 10.0);

# ╔═╡ bb287258-f383-4c04-aee2-c86b8a39a9bd
plot(medium2D, fields=[:vp, :rho])

# ╔═╡ Cell order:
# ╟─86d3f068-a979-42f5-a9e7-138e94c16b38
# ╟─de2f97fa-f64d-4754-bd34-e44dbf13c336
# ╟─0ffd8ee4-1735-4756-befb-c7c10d08eb34
# ╟─7687d367-f9d5-4539-9156-e26d87379f87
# ╟─d9b71485-8b64-4ad0-a242-dde6300af835
# ╠═981f55af-1557-49c5-921d-2e7e343a511b
# ╠═a6b59ba1-ea71-467f-bf26-49a634a1bbd9
# ╠═3c549ecd-52dc-47d2-b846-620148089296
# ╟─addd0bef-5d52-4e06-9350-f6dade3db428
# ╟─50494fe6-6bdf-4db1-8e61-6c9b39627b1f
# ╠═80109fc0-7446-47eb-b912-1f767d824fbd
# ╠═701d36a4-5f09-43ab-bd09-aefc92196182
# ╠═63ec2c93-539e-4a58-a10a-fb61602ae03e
# ╠═d356c66b-b07e-415a-a763-9d7f235f1c48
# ╠═bdac654f-5743-477c-af43-17cb8ffc9b60
# ╠═3cc49ee4-06f6-4459-9f42-8838963f344a
# ╠═15b42a58-5ad7-42a6-9dd4-7db648c33f55
# ╠═fc58c853-2665-45f8-997d-1b3217fa8ec0
# ╠═7a721ee1-227b-4975-b4c7-1ac1ec83e973
# ╠═e55ed529-8670-4326-ac26-e3eaa55f624b
# ╠═566a79e2-2e79-40db-8a46-e0a341989bb0
# ╟─bc26144a-167b-4e99-a631-d7e16381fe97
# ╠═3b9d77ed-ad38-429d-b6de-34722a3850a8
# ╠═5bf36510-ddb5-43bb-a3a0-02a22578b27a
# ╠═315f5412-74ec-4b18-9286-9201c7ab5dc1
# ╠═bb287258-f383-4c04-aee2-c86b8a39a9bd
