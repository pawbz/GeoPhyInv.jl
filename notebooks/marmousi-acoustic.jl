### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ 28f0f643-aa80-4f92-a0f1-10cf88be3383
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

# ╔═╡ ca53ee61-98ad-436a-bd10-14290c0ca01e
begin
    Pkg.activate(Base.current_project())
    Pkg.instantiate()
end

# ╔═╡ 58ec85a4-792a-4294-9147-e3c9637e466e
@revise using GeoPhyInv

# ╔═╡ 2520f1ee-a102-4b08-9ebd-1b9668346569
TableOfContents()

# ╔═╡ df39c4c2-3be8-4675-ad8d-1262c280bbc3
Pkg.activate

# ╔═╡ 6d1c8a37-b872-4041-8544-a3437ebe01a1
marm=ElasticMedium(Marmousi2())

# ╔═╡ b4ba72d6-0916-4738-8d25-c807e636ce2b
ageom = AGeom(marm.grid, :surf, SSrcs(1), Recs(100));

# ╔═╡ 0c8c94c7-c2aa-457f-9ab2-2f80935b1534
plot(marm, ageom, fields=[:vp, :rho])

# ╔═╡ 45bf3357-3e14-483f-acc4-ed3df12fa9fb
wav, tgrid = ricker(marm, 10, 2)

# ╔═╡ 2095e9a7-a6bc-4399-989b-4c57b3de7afd
srcwav = Srcs(tgrid, ageom, [:p])

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
# ╠═2520f1ee-a102-4b08-9ebd-1b9668346569
# ╠═28f0f643-aa80-4f92-a0f1-10cf88be3383
# ╠═df39c4c2-3be8-4675-ad8d-1262c280bbc3
# ╠═ca53ee61-98ad-436a-bd10-14290c0ca01e
# ╠═58ec85a4-792a-4294-9147-e3c9637e466e
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
