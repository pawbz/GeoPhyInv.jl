### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ 40eedc0e-f934-4dd4-b00c-34d47a027dd8
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

# ╔═╡ 6520bf6f-332b-4712-8fec-1d824114c4ef
begin
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
end

# ╔═╡ afbd2318-5edf-46a1-83ad-49ed37247810
@revise using GeoPhyInv

# ╔═╡ b48a5c89-775c-426a-926d-2c6136112d31
TableOfContents()

# ╔═╡ 24e8c3e4-c42d-4fab-9a57-875fcae1c4bc
md"# Gallery"

# ╔═╡ 179140e9-9175-44e1-80c5-5c5b4d33ac9e
md"## Medium"

# ╔═╡ 778490ec-c97b-4587-9508-38c1fdc6539e
medium1 = AcousticMedium(Homogeneous(), 5.0)

# ╔═╡ cc60f80a-9de8-4568-80a6-f8e6a8fdee97
medium2 = ElasticMedium(Homogeneous(), 5.0)

# ╔═╡ eb5fcc15-c076-4d0f-8f53-b5ba4d826063
medium1.mgrid

# ╔═╡ d77093ed-e18b-4d4d-b3b4-2011fe0fa465
plot(medium1)

# ╔═╡ d20aaa97-e5f5-43d6-8af4-6050c7f1a48b
plot(medium2)

# ╔═╡ 00225323-b954-4a80-a1bb-4306bc899124
typeof(medium) |> subtype

# ╔═╡ 0affb9ec-7607-4d87-a2e9-0f4b03411b4b
a=fill(Float32(2500.0), 10,10)

# ╔═╡ 0fda6f68-3ae4-4e26-b7d5-3136760ebd57
convert(eltype(a), 5.0)

# ╔═╡ d8ab1383-0aed-417a-b6b1-a25a4db67990
Homogeneous()

# ╔═╡ 4e442a17-028e-4612-a8e6-0e22811813e5
Marmousi2()

# ╔═╡ c908bae1-9ee4-4684-a53b-a83e8137201a
fill(Float32(4.0), 4,5)

# ╔═╡ d2ec4b5e-0368-4e31-a2d7-73e79d576e95
ndims(medium2)

# ╔═╡ 2d97c59c-69a3-47e2-bf1c-67bd2abe2ce3
md"## SeisForwExpt"

# ╔═╡ 5fb98ed2-3c57-45d8-a610-55d532e9235c
SeisForwExpt(FdtdAcoustic(), Homogeneous())

# ╔═╡ 60f3f582-9e47-40d6-8ecd-a3a6603a6394
SeisForwExpt(FdtdAcoustic(), RandScatterer())

# ╔═╡ 9fa2ff76-774a-4d8e-9d1a-51c54cb39797


# ╔═╡ Cell order:
# ╠═b48a5c89-775c-426a-926d-2c6136112d31
# ╠═40eedc0e-f934-4dd4-b00c-34d47a027dd8
# ╠═6520bf6f-332b-4712-8fec-1d824114c4ef
# ╠═afbd2318-5edf-46a1-83ad-49ed37247810
# ╠═24e8c3e4-c42d-4fab-9a57-875fcae1c4bc
# ╟─179140e9-9175-44e1-80c5-5c5b4d33ac9e
# ╠═778490ec-c97b-4587-9508-38c1fdc6539e
# ╠═cc60f80a-9de8-4568-80a6-f8e6a8fdee97
# ╠═eb5fcc15-c076-4d0f-8f53-b5ba4d826063
# ╠═d77093ed-e18b-4d4d-b3b4-2011fe0fa465
# ╠═d20aaa97-e5f5-43d6-8af4-6050c7f1a48b
# ╠═00225323-b954-4a80-a1bb-4306bc899124
# ╠═0affb9ec-7607-4d87-a2e9-0f4b03411b4b
# ╠═0fda6f68-3ae4-4e26-b7d5-3136760ebd57
# ╠═d8ab1383-0aed-417a-b6b1-a25a4db67990
# ╠═4e442a17-028e-4612-a8e6-0e22811813e5
# ╠═c908bae1-9ee4-4684-a53b-a83e8137201a
# ╠═d2ec4b5e-0368-4e31-a2d7-73e79d576e95
# ╠═2d97c59c-69a3-47e2-bf1c-67bd2abe2ce3
# ╠═5fb98ed2-3c57-45d8-a610-55d532e9235c
# ╠═60f3f582-9e47-40d6-8ecd-a3a6603a6394
# ╠═9fa2ff76-774a-4d8e-9d1a-51c54cb39797
