### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ e9e81c62-bef9-4736-9fce-bfe3b0cada7a
using SparseArrays, Test

# ╔═╡ c23c557c-a426-48d0-a1dd-7cf2f3a47fa5
"""
minimum index using imask
"""
function indminimum(x, val, imask=[])
    min_i = 0
    min_x = typemax(Float64)
    for (i, xi) in enumerate(x)
        dist = abs(xi - val)
        if ((dist < min_x) && (i ∉ imask))
            min_x = dist
            min_i = i
        end
    end
    return min_i
end


# ╔═╡ 8d41c81f-390d-417e-9231-7b42184dfb18
"""
Return n indices in order
Cannot find a julia method which does, this.
If a faster method is found, replace it later.
"""
function indminn!(ivec, x, val)
    # using enumerate to avoid indexing
    n = length(ivec)
    if (length(x) < n)
        # if the desired indices is greater than length of x, what to do..
        fill!(ivec, one(eltype(ivec)))
    else
        fill!(ivec, zero(eltype(ivec)))
        for inn in 1:n
            ivec[inn] = indminimum(x, val, ivec)
        end
        sort!(ivec)
    end
    return ivec
end


# ╔═╡ 6d953eda-506e-484c-b084-5c318c8f67d0
function indminn(x, val, n)
    ivec = fill(0, n)
    indminn!(ivec, x, val)
    return ivec
end

# ╔═╡ 7dd272c1-f614-4410-9aba-192e2a33a85a
"""
Output CartesianIndices, LinearIndices and corresponding weights for linear interpolation of point P on mgrid.
mgrid is created depending on the input field.
For example, when 2D,
mgrid=[range(1.3, stop=10.6,step=0.003), range(1.2,stop=15.3,step=0.004)]
P=[5,5]

I, J, L, V are used to create sparse spray and interpolation matrices for sources and receivers.
"""
function find_neighbour_weights(mgrid, P)
    @assert length(mgrid) == length(P)
    N = length(mgrid)
    idx = [indminn(mgrid[i], P[i], 2) for i = 1:N]
    denomI = inv(prod([diff(mgrid[i][idx[i]])[1] for i = 1:N]))
    c = CartesianIndices(Tuple(broadcast(x -> x[1]:x[2], idx)))
    l = LinearIndices(Tuple(length.(mgrid)))
    # linear indices 
    L = [l[cc] for cc in c]
    # calculate distance to each neighbour
    dist = [sum([abs2(mgrid[i][cc[i]] - P[i]) for i = 1:N]) * denomI for cc in c]
    # the weight is inversely proportional to the distance
    Idist = sortperm(vec(dist), rev=true)
    # the weights are proportional to the area patches
    V = [prod([abs(mgrid[i][cc[i]] - P[i]) for i = 1:N]) * denomI for cc in c]
    return vec(L), vec(V)[Idist]
end

# ╔═╡ 1c90e052-ed95-11ed-0cad-ff4e7234f520
function get_proj_matrix(migrid::T, mmgrid::T; use_gpu=_fd_use_gpu, number=Data.Number) where {T}
    ijv = mapreduce((x, y) -> (; I=vcat(x.I, y.I), J=vcat(x.J, y.J), V=vcat(x.V, y.V)), 1:prod(length.(migrid)), Iterators.product(migrid...)) do iP, P
        Ln, Vn = find_neighbour_weights(mmgrid, P)
        return (; I=Ln, J=fill(iP, length(Ln)), V=number.(Vn))
    end
    (; I, J, V) = ijv

    if (use_gpu)
        ssprayw = CuSparseMatrixCSC(sparse(I, J, V, prod(length.(mmgrid)), prod(length.(migrid))))
    else
        ssprayw = sparse(I, J, V, prod(length.(mmgrid)), prod(length.(migrid)))
    end
end

# ╔═╡ 137e4dfd-60f0-460b-8239-082680b09b77
function get_proj_matrix(Ps, mmgrid; use_gpu=_fd_use_gpu, number=Data.Number)
    ijv = mapreduce((x, y) -> (; I=vcat(x.I, y.I), J=vcat(x.J, y.J), V=vcat(x.V, y.V)), 1:length(Ps), Ps) do iP, P
        Ln, Vn = find_neighbour_weights(mmgrid, P)
        return (; I=Ln, J=fill(iP, length(Ln)), V=number.(Vn))
    end
    (; I, J, V) = ijv

    if (use_gpu)
        ssprayw = CuSparseMatrixCSC(sparse(I, J, V, prod(length.(mmgrid)), length(Ps)))
    else
        ssprayw = sparse(I, J, V, prod(length.(mmgrid)), length(Ps))
    end
end

# ╔═╡ 205e6872-a5fa-4b2d-9d1d-8f0c55793f30
get_proj_matrix( [[1.2], [1.5]],[1:10], use_gpu=false, number=Float32)

# ╔═╡ d0de00fa-9fa7-4159-9756-97e4c4bc33e2
@time get_proj_matrix([range(-0.75, 0.5, length=3)], [range(-1, 1, length=10)], use_gpu=false, number=Float32)'

# ╔═╡ fdef69e3-41cf-4bf7-b50e-3ed097c58818
get_proj_matrix([range(-1, 1, length=10)], [range(-0.75, 0.5, length=3)], use_gpu=false, number=Float32)

# ╔═╡ eecdb64c-c123-4678-b449-66cdef070244
get_proj_matrix([range(-0.75, 0.5, length=3), range(-0.75, 0.5, length=3)], [range(-1, 1, length=10), range(-1, 1, length=10)], use_gpu=false, number=Float32)'

# ╔═╡ 0c7520d9-60a9-4925-8b58-354ea863ecca
find_neighbour_weights([1:10], [1.2])

# ╔═╡ 95ff1284-db71-4927-b509-e4f76da8c892
@test sum(find_neighbour_weights([1:10, 1:10], [2.2, 4.2])[2]) == 1

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "4d18121a8e51d2800755224622a32a12acdf6109"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═e9e81c62-bef9-4736-9fce-bfe3b0cada7a
# ╠═1c90e052-ed95-11ed-0cad-ff4e7234f520
# ╠═137e4dfd-60f0-460b-8239-082680b09b77
# ╠═7dd272c1-f614-4410-9aba-192e2a33a85a
# ╠═0c7520d9-60a9-4925-8b58-354ea863ecca
# ╠═205e6872-a5fa-4b2d-9d1d-8f0c55793f30
# ╠═95ff1284-db71-4927-b509-e4f76da8c892
# ╠═d0de00fa-9fa7-4159-9756-97e4c4bc33e2
# ╠═fdef69e3-41cf-4bf7-b50e-3ed097c58818
# ╠═eecdb64c-c123-4678-b449-66cdef070244
# ╠═c23c557c-a426-48d0-a1dd-7cf2f3a47fa5
# ╠═6d953eda-506e-484c-b084-5c318c8f67d0
# ╠═8d41c81f-390d-417e-9231-7b42184dfb18
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
