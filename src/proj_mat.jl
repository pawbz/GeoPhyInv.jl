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
function find_neighbour_weights(mgrid, P, weight_fn=(w, cc) -> w)
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
    V = [weight_fn(prod([abs(mgrid[i][cc[i]] - P[i]) for i = 1:N]) * denomI, cc) for cc in c]
    return vec(L), vec(V)[Idist]
end

# ╔═╡ 1c90e052-ed95-11ed-0cad-ff4e7234f520
function get_proj_matrix(migrid, mmgrid; use_gpu=_fd_use_gpu, weight_fn=(w, cc) -> w)
    ijv = mapreduce((x,y)->(; I=vcat(x.I,y.I), J=vcat(x.J,y.J), V=vcat(x.V, y.V)), 1:prod(length.(migrid)), Iterators.product(migrid...)) do iP, P
        Ln, Vn = find_neighbour_weights(mmgrid, P, weight_fn)
        return (; I=Ln, J=fill(iP, length(Ln)), V=Vn)
    end
	(; I, J, V) = ijv

    if (use_gpu)
        ssprayw = CuSparseMatrixCSC(sparse(I, J, V, prod(length.(mmgrid)), prod(length.(migrid))))
    else
        ssprayw = sparse(I, J, V, prod(length.(mmgrid)), prod(length.(migrid)))
    end
end

# ╔═╡ d0de00fa-9fa7-4159-9756-97e4c4bc33e2
@time get_proj_matrix([range(-0.75, 0.5, length=3)], [range(-1, 1, length=10)], use_gpu=false)'

# ╔═╡ eecdb64c-c123-4678-b449-66cdef070244
get_proj_matrix([range(-0.75, 0.5, length=3), range(-0.75, 0.5, length=3)], [range(-1, 1, length=10), range(-1, 1, length=10)], use_gpu=false)'

# ╔═╡ fdef69e3-41cf-4bf7-b50e-3ed097c58818
get_proj_matrix([range(-1, 1, length=10)], [range(-0.75, 0.5, length=3)], use_gpu=false)

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

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "4d18121a8e51d2800755224622a32a12acdf6109"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"
"""

# ╔═╡ Cell order:
# ╠═e9e81c62-bef9-4736-9fce-bfe3b0cada7a
# ╠═1c90e052-ed95-11ed-0cad-ff4e7234f520
# ╠═7dd272c1-f614-4410-9aba-192e2a33a85a
# ╠═0c7520d9-60a9-4925-8b58-354ea863ecca
# ╠═95ff1284-db71-4927-b509-e4f76da8c892
# ╠═d0de00fa-9fa7-4159-9756-97e4c4bc33e2
# ╠═eecdb64c-c123-4678-b449-66cdef070244
# ╠═fdef69e3-41cf-4bf7-b50e-3ed097c58818
# ╠═c23c557c-a426-48d0-a1dd-7cf2f3a47fa5
# ╠═6d953eda-506e-484c-b084-5c318c8f67d0
# ╠═8d41c81f-390d-417e-9231-7b42184dfb18
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
