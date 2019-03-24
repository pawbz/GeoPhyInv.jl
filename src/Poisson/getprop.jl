"""
Return stuff stored in `PoissonExpt` after using `mod!`. 
"""
function Base.getindex(pa::ParamExpt, s::Symbol)
	@assert s in [:data]
	if(s==:data)
		return pa.data
	end
end


