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

# ╔═╡ 0ffd8ee4-1735-4756-befb-c7c10d08eb34
# ╠═╡ show_logs = false
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

# ╔═╡ 0285e307-c619-4fa1-9e12-0ba755d476fb
begin
    Pkg.activate(Base.current_project())
    Pkg.instantiate()
end

# ╔═╡ 59743052-ca23-4a0d-a32d-e66adb6cd352
@revise using GeoPhyInv

# ╔═╡ d9b71485-8b64-4ad0-a242-dde6300af835
# ╠═╡ show_logs = false

using Statistics, LossFunctions, LinearAlgebra


# ╔═╡ 45656ece-32e9-488f-be20-6546017e1e94
TableOfContents()

# ╔═╡ 2c5e29cf-67cd-4c74-b7e8-8f16b1828390
pa_acoustic = SeisForwExpt(:acou_homo2D);

# ╔═╡ 29ad7daf-873f-4b59-b10b-780d6798426e


# ╔═╡ b236913b-5320-461a-8583-74eb4140ff27
pa_elastic = SeisForwExpt(:elastic_homo2D);

# ╔═╡ 3109122c-5a01-43b7-9ade-46deeff185a7
@bind src_type Radio(["1" => 1, "2" => 2], default="1")

# ╔═╡ 309e105e-74ca-427a-aca3-bc7212eb88fa
@bind src_field Radio(["p" => "p", "vz" => "vz", "vx" => "vx"], default="vz")

# ╔═╡ 4a153cba-fb5c-40d8-a82f-1dbbc1d7f6c6
src_type

# ╔═╡ bf23dad7-2c0a-47f1-ad52-7f38a3658ae1
function test_backprop(pa, fields)
    println("############ Testing Backprop for source type ", src_field)

    srcwav = pa[:srcwav][1]
    GeoPhyInv.update!(srcwav[1], [Symbol(src_field)])

    # save boundary and final states
    pa.c.attrib_mod.mode = :forward

    GeoPhyInv.update!(pa, [srcwav], [parse(Int, src_type)])

    update!(pa)

	snaps_forw = deepcopy(pa[:snaps, 1])


    # change source flag and update wavelets in pa
    GeoPhyInv.update!(pa, [srcwav], [-parse(Int, src_type)])

    # force boundary values and use saved initial state
    pa.c.attrib_mod.mode = :adjoint

    update!(pa)
		
	snaps_back = deepcopy(pa[:snaps, 1])
  return snaps_forw, snaps_back

end

# ╔═╡ c4e9f408-8408-41b7-a788-20e50dace617
snaps_forw, snaps_back = test_backprop(pa_acoustic, [:p, :vx, :vz])


# ╔═╡ 4bd8806f-5afe-48af-83cd-aa8f6220eb44
compute L2dist
    @show err = value(L2DistLoss(), vec(rec1), vec(rec2), AggMode.Mean())

    # desired accuracy?
    @test err < 1e-10

# ╔═╡ Cell order:
# ╠═45656ece-32e9-488f-be20-6546017e1e94
# ╠═0ffd8ee4-1735-4756-befb-c7c10d08eb34
# ╠═0285e307-c619-4fa1-9e12-0ba755d476fb
# ╠═59743052-ca23-4a0d-a32d-e66adb6cd352
# ╠═d9b71485-8b64-4ad0-a242-dde6300af835
# ╠═2c5e29cf-67cd-4c74-b7e8-8f16b1828390
# ╠═c4e9f408-8408-41b7-a788-20e50dace617
# ╠═29ad7daf-873f-4b59-b10b-780d6798426e
# ╠═b236913b-5320-461a-8583-74eb4140ff27
# ╠═3109122c-5a01-43b7-9ade-46deeff185a7
# ╠═309e105e-74ca-427a-aca3-bc7212eb88fa
# ╠═4a153cba-fb5c-40d8-a82f-1dbbc1d7f6c6
# ╠═bf23dad7-2c0a-47f1-ad52-7f38a3658ae1
# ╠═4bd8806f-5afe-48af-83cd-aa8f6220eb44
