"""
Return snaps stored in `SeisForwExpt` after using `mod!`. 

# Arguments

* `iss::Int64=1` : supersource index
"""
function Base.getindex(pa::Param, s::Symbol, iss::Int=1)
	@assert s in [:snaps]
	if(s==:snaps)
		return pa.p.localpart.ss[iss].snaps
	end
end


