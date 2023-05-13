### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ e9e81c62-bef9-4736-9fce-bfe3b0cada7a
using SparseArrays, Test

# ╔═╡ 78c2e09e-843c-4b32-a09a-df4027d10895
function get_neighbour_indices(arr, val)
    len = length(arr)
    idx = searchsortedfirst(arr, val)
    if idx == 1
        return [1, 2]
	elseif idx >= len
        return [len-1, len]
    else
        return [idx-1, idx]
    end
end

# ╔═╡ 49d5dc87-0256-404c-9b59-f1db4450a3e3
"""
This function takes in three vectors x, y, z, and three values xi, yi, zi, corresponding to the interpolation point. It first checks if the interpolation point is inside the grid, then finds the indices of the 8 neighbouring grid points using the get_neighbour_indices function. It then computes the interpolation weights and stores them in a sparse vector with values of 1 and the specified weights at the corresponding indices. The resulting sparse vector has length n * m * p and a single column.
"""
function bilinear_interp(x, y, z, xi, yi, zi; number=Data.Number)
    if (xi < minimum(x) || xi > maximum(x) || yi < minimum(y) || yi > maximum(y) || zi < minimum(z) || zi > maximum(z))
        throw(ArgumentError("Interpolation point is outside the grid."))
    end
    n, m, p = length(x), length(y), length(z)
    x_idx1, x_idx2 = get_neighbour_indices(x, xi)
    y_idx1, y_idx2 = get_neighbour_indices(y, yi)
    z_idx1, z_idx2 = get_neighbour_indices(z, zi)
    l = LinearIndices((n, m, p))
    i000 = l[x_idx1, y_idx1, z_idx1]
    i001 = l[x_idx1, y_idx1, z_idx2]
    i010 = l[x_idx1, y_idx2, z_idx1]
    i011 = l[x_idx1, y_idx2, z_idx2]
    i100 = l[x_idx2, y_idx1, z_idx1]
    i101 = l[x_idx2, y_idx1, z_idx2]
    i110 = l[x_idx2, y_idx2, z_idx1]
    i111 = l[x_idx2, y_idx2, z_idx2]
    dx = (xi - x[x_idx1]) / (x[x_idx2] - x[x_idx1])
    dy = (yi - y[y_idx1]) / (y[y_idx2] - y[y_idx1])
    dz = (zi - z[z_idx1]) / (z[z_idx2] - z[z_idx1])
    w000 = (1 - dx) * (1 - dy) * (1 - dz)
    w001 = (1 - dx) * (1 - dy) * dz
    w010 = (1 - dx) * dy * (1 - dz)
    w011 = (1 - dx) * dy * dz
    w100 = dx * (1 - dy) * (1 - dz)
    w101 = dx * (1 - dy) * dz
    w110 = dx * dy * (1 - dz)
    w111 = dx * dy * dz
    weights = [w000, w001, w010, w011, w100, w101, w110, w111]
    indices = [i000, i001, i010, i011, i100, i101, i110, i111]
    return sparse(indices, fill(1, 8), number.(weights), n * m * p, 1)
end

# ╔═╡ f8a70920-ac87-4126-980f-c3cbd223269a
"""
This function takes as input the x and y coordinates of the grid points, the corresponding function values z, and the x and y coordinates of the point to interpolate (xi, yi). It returns a sparse matrix with the weights for each of the four grid points used in the bilinear interpolation.

Note that if the interpolation point is outside the grid, the function will throw an error.
"""
function bilinear_interp(x, y, xi, yi; number=Data.Number)
	if (xi < minimum(x) || xi > maximum(x) || yi < minimum(y) || yi > maximum(y))
        throw(ArgumentError("Interpolation point is outside the grid."))
    end
    n, m = length(x), length(y)
    x_idx1, x_idx2 = get_neighbour_indices(x, xi)
	y_idx1, y_idx2 = get_neighbour_indices(y, yi)
	
    l = LinearIndices((n, m))
	
    i00 = l[x_idx1, y_idx1]
    i01 = l[x_idx1, y_idx2]
    i10 = l[x_idx2, y_idx1]
    i11 = l[x_idx2, y_idx2]
    
    dx = (xi - x[x_idx1]) / (x[x_idx2] - x[x_idx1])
    dy = (yi - y[y_idx1]) / (y[y_idx2] - y[y_idx1])
    
    w00 = (1-dx)*(1-dy)
    w01 = (1-dx)*dy
    w10 = dx*(1-dy)
    w11 = dx*dy
    
    weights = [w00, w01, w10, w11]
    indices = [i00, i01, i10, i11]
    
    return sparse(indices, fill(1, 4), number.(weights), n*m, 1)
end

# ╔═╡ 4fa99e86-d6fb-4da1-91d9-8065245c1191
"""
This function returns a sparse matrix with four nonzero elements, which are the weights for the four surrounding points used in bilinear interpolation. If the given xi value is outside the range of x, the function returns a sparse matrix with a single nonzero element, which is the index of the closest boundary point.
"""
function bilinear_interp(x, xi; number=Data.Number)
	if (xi < minimum(x) || xi > maximum(x))
        throw(ArgumentError("Interpolation point is outside the grid."))
    end
    n = length(x)
    x_idx1, x_idx2 = get_neighbour_indices(x, xi)

    l = LinearIndices((n,))
	
    i00 = l[x_idx1]
    i01 = l[x_idx2]
   
    dx = (xi - x[x_idx1]) / (x[x_idx2] - x[x_idx1])
    
    w00 = 1 - dx
    w01 = dx
    
    weights = [w00, w01]
    indices = [i00, i01]
    
    return sparse(indices, fill(1, 2), number.(weights), n, 1)
end

# ╔═╡ 1c90e052-ed95-11ed-0cad-ff4e7234f520
function get_proj_matrix(migrid::T, mmgrid::T; use_gpu=_fd_use_gpu, number=Data.Number) where {T}
    mat = mapreduce(sparse_hcat, Iterators.product(migrid...)) do P
  		bilinear_interp(mmgrid..., P..., number=number)
    end
    if (use_gpu)
        return CuSparseMatrixCSC(mat)
    else
        return mat
    end
end

# ╔═╡ 137e4dfd-60f0-460b-8239-082680b09b77
function get_proj_matrix(Ps, mmgrid; use_gpu=_fd_use_gpu, number=Data.Number)
    mat = mapreduce(sparse_hcat, Ps) do P
		bilinear_interp(mmgrid..., P..., number=number)
    end
    if (use_gpu)
        return CuSparseMatrixCSC(mat)
    else
        return mat
    end
end

# ╔═╡ 205e6872-a5fa-4b2d-9d1d-8f0c55793f30
get_proj_matrix( [[1.2], [1.5]],[1:10], use_gpu=false, number=Float32)

# ╔═╡ d0de00fa-9fa7-4159-9756-97e4c4bc33e2
@time get_proj_matrix([range(-0.75, 0.5, length=3)], [range(-1, 1, length=10)], use_gpu=false, number=Float32)

# ╔═╡ eecdb64c-c123-4678-b449-66cdef070244
get_proj_matrix([range(-0.75, 0.5, length=3), range(-0.75, 0.5, length=3)], [range(-1, 1, length=10), range(-1, 1, length=10)], use_gpu=false, number=Float32)'

# ╔═╡ e2416fd8-e0be-49db-a672-0d9fd9b09db4
@time bilinear_interp(range(-1, 1, length=5), range(-1, 1, length=5), 0.5, 0.2)

# ╔═╡ be183231-697e-408f-aa08-ad0d6db1f663
@time bilinear_interp(1:10, 1.8, number=Float32)

# ╔═╡ 95ff1284-db71-4927-b509-e4f76da8c892
@test sum(bilinear_interp(1:10, range(10, 20, length=7), 2.2, 14.4, number=Float64)) ≈ 1

# ╔═╡ e1c0c307-2fb4-4fe8-a6d0-cb26624c9075
@test sum(bilinear_interp(1:10, range(10, 20, length=7), range(30, 41, length=5), 5.2, 14.2, 35.0, number=Float64)) ≈ 1

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
# ╠═49d5dc87-0256-404c-9b59-f1db4450a3e3
# ╠═f8a70920-ac87-4126-980f-c3cbd223269a
# ╠═e2416fd8-e0be-49db-a672-0d9fd9b09db4
# ╠═4fa99e86-d6fb-4da1-91d9-8065245c1191
# ╠═be183231-697e-408f-aa08-ad0d6db1f663
# ╠═205e6872-a5fa-4b2d-9d1d-8f0c55793f30
# ╠═95ff1284-db71-4927-b509-e4f76da8c892
# ╠═e1c0c307-2fb4-4fe8-a6d0-cb26624c9075
# ╠═d0de00fa-9fa7-4159-9756-97e4c4bc33e2
# ╠═eecdb64c-c123-4678-b449-66cdef070244
# ╠═78c2e09e-843c-4b32-a09a-df4027d10895
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
