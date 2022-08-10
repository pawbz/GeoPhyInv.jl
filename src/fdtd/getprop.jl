"""
An instance of `SeisForwExpt`. 

## Indexing
* `pa[:data, i]` : get data for ith supersource (defaults to first), after `update!`
* `pa[:ageom]` : get `ageom` originally input while creating `pa`
* `pa[:srcwav]`: `srcwav` used to create `pa`
* `pa[:snaps, i]` : snapshots corresponding to `tsnaps` and ith source, after `update!`
* `pa[:medium]`: medium 
* `pa[:exmedium]`: extended medium to PML regions
"""
function Base.getindex(pa::PFdtd, s::Symbol, iss::Int = 1)
    @assert iss ≤ length(pa.c.ageom[1])
    @assert s in vcat([:snaps, :data], collect(fieldnames(P_common)))
    if (s == :snaps)
        ip = findall(in.(iss, pa.sschunks))[1]
        issp = findall(x -> x == iss, pa.sschunks[ip])[]
        @sync begin
            p = procs(pa.p)[ip]
            @sync remotecall_wait(p) do
                pap = localpart(pa.p)
                return pap[1].ss[issp].snaps
            end
        end
    elseif (s == :data)
        (iss == 0) ? (return pa.c.data[1]) : (return pa.c.data[1][iss])
    else
        return getfield(pa.c, s)
    end
end

function Base.show(io::Base.IO, pa::PFdtd)
    return println(@doc Base.getindex(::GeoPhyInv.PFdtd, ::Symbol, ::Int))
end

Base.show(io::Base.IO, pa::Vector{P_x_worker_x_pw{N}}) where {N} = nothing
Base.show(io::Base.IO, pa::P_x_worker_x_pw) = nothing
Base.show(io::Base.IO, pa::P_x_worker_x_pw_x_ss) = nothing



