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

# ╔═╡ 19ca0fd6-a97e-47e2-a950-d7e8a6ba8b88
begin
    import Pkg
    Pkg.add("PlutoLinks")
    Pkg.add("PlutoUI")
    Pkg.add("PlutoTest")
    Pkg.add("Plots")
    using PlutoLinks: @revise
    using PlutoUI, PlutoTest, Plots
end

# ╔═╡ 00333dab-a753-4073-ae6a-36aee85b93c9
TableOfContents()

# ╔═╡ 7ceb8cb6-5976-4efd-aed9-d54772968de2
@bind reload_geophyinv Button("using GeoPhyInv")

# ╔═╡ e2ddc23f-99f8-405e-8846-dd074e306630
begin
    reload_geophyinv
    Pkg.activate(Base.current_project())
    Pkg.instantiate()
    @revise using GeoPhyInv
end

# ╔═╡ d61191f0-4bc4-466b-995c-96a137056ccb
begin
    reload_geophyinv
    GeoPhyInv.@init_parallel_stencil(2, false, Float32, 2)
	using Statistics
end

# ╔═╡ 3b47fd1e-3d24-4346-884d-4e0c0db228cf
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

# ╔═╡ 6d1c8a37-b872-4041-8544-a3437ebe01a1
marm=Medium(:marmousi2)

# ╔═╡ b4ba72d6-0916-4738-8d25-c807e636ce2b
ageom = AGeom(marm.mgrid, :surf, SSrcs(1), Recs(100));

# ╔═╡ 0c8c94c7-c2aa-457f-9ab2-2f80935b1534
plot(marm, ageom, fields=[:vp, :rho])

# ╔═╡ 45bf3357-3e14-483f-acc4-ed3df12fa9fb
wav, tgrid = ricker(marm, 10, 2)

# ╔═╡ 2095e9a7-a6bc-4399-989b-4c57b3de7afd
srcwav = SrcWav(tgrid, ageom, [:p])

# ╔═╡ ee9aa323-80f0-4bf3-b82d-df572cf58174
update!(srcwav, [:p], wav)

# ╔═╡ bb69a19c-aef2-4fa7-a314-d5f00cd6c38f
plot(srcwav[1])

# ╔═╡ 94683a46-8755-473c-b4f1-3fe616f44012
pa_acoustic = SeisForwExpt(
    FdtdAcoustic(),
    medium = marm,
    ageom = ageom,
	tgrid = tgrid,
    srcwav = srcwav,
	rfields = [:p],
   );

# ╔═╡ 1f0a5b21-0971-421e-8fe9-b6e658cb090d
update!(pa_acoustic)

# ╔═╡ be44bca9-88fa-4549-9f89-6391cef1d9f7
plot(pa_acoustic[:data], 99.9)

# ╔═╡ Cell order:
# ╟─00333dab-a753-4073-ae6a-36aee85b93c9
# ╟─7ceb8cb6-5976-4efd-aed9-d54772968de2
# ╟─3b47fd1e-3d24-4346-884d-4e0c0db228cf
# ╟─19ca0fd6-a97e-47e2-a950-d7e8a6ba8b88
# ╟─e2ddc23f-99f8-405e-8846-dd074e306630
# ╠═d61191f0-4bc4-466b-995c-96a137056ccb
# ╠═6d1c8a37-b872-4041-8544-a3437ebe01a1
# ╠═b4ba72d6-0916-4738-8d25-c807e636ce2b
# ╠═0c8c94c7-c2aa-457f-9ab2-2f80935b1534
# ╠═45bf3357-3e14-483f-acc4-ed3df12fa9fb
# ╠═2095e9a7-a6bc-4399-989b-4c57b3de7afd
# ╠═ee9aa323-80f0-4bf3-b82d-df572cf58174
# ╠═bb69a19c-aef2-4fa7-a314-d5f00cd6c38f
# ╠═94683a46-8755-473c-b4f1-3fe616f44012
# ╠═1f0a5b21-0971-421e-8fe9-b6e658cb090d
# ╠═be44bca9-88fa-4549-9f89-6391cef1d9f7
