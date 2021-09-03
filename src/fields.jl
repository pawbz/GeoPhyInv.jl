# Stress and velocity fields (used for multiple dispatch)


"""
Return axis names of 1D, 2D or 3D fields
"""
function dim_names(N,prefix="",suffix="";order=1,)
	if(N==3)
		names=(order==2) ? [:zz,:yy,:xx,:xz,:xy,:yz] : [:z,:y,:x]
	elseif(N==2)
		names=(order==2) ? [:zz,:xx,:xz] : [:z,:x]
	elseif(N==1)
		names=[:z]
	else
		error("invalid dim num")
	end
	return [Symbol(prefix, string(m), suffix) for m in names]
end


# pressure for acoustic 
struct p end 
# stress and velocity (tauxx, tauxy, tauxz,)
for f in vcat(dim_names(3,"v"), dim_names(3,"tau";order=2))
	@eval struct $f end
end

# tauii sam
function zeros(::$f,nz,ny,nx)
	
end