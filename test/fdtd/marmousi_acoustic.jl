### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 87ef6a22-0fda-11ed-37f7-c129450c0930
import Pkg; Pkg.activate("GeoPhyInv")

# ╔═╡ 76b7a72a-0fdc-47c3-9ee0-3fc20e86526a
begin
	using Revise
	using GeoPhyInv
	using Statistics
	using Test
	using CUDA
	using PlutoUI
	CUDA.device!(1)
end

# ╔═╡ eb89564f-b1c5-47e8-9899-fa4bf22abc64
using Plots; gr()

# ╔═╡ c8118bea-1fc1-4d67-9ea6-54a5e6127866
@init_parallel_stencil(2, true, Float32, 2)

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
# ╠═87ef6a22-0fda-11ed-37f7-c129450c0930
# ╠═76b7a72a-0fdc-47c3-9ee0-3fc20e86526a
# ╠═eb89564f-b1c5-47e8-9899-fa4bf22abc64
# ╠═c8118bea-1fc1-4d67-9ea6-54a5e6127866
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
