
"""
## Indexing
`pa` is an instance of `SeisForwExpt`. 
* `pa[:data,i]` : get data for ith supersource (defaults to first), after `update!`
* `pa[:ageom]` : get `ageom` originally input while creating `pa`
* `pa[:srcwav]`: `srcwav` used to create `pa`
"""
function Base.getindex(pa::PBorn, s::Symbol, iss::Int = 1)
    @assert iss ≤ length(pa.ageom)
    @assert s in [:data, :ageom, :srcwav]
    if (s == :data)
        (iss == 0) ? (return pa.data) : (return pa.data[iss])
    elseif (s == :ageom)
        return pa.ageom
    elseif (s == :srcwav)
        return pa.srcwav
    end
end

function Base.show(io::Base.IO, pa::PFdtd)
    println(typeof(pa), "")
    println("pa[:data] : modeled data after running `update!`")
    println("pa[:ageom]")
    println("pa[:srcwav]")
end
