

"""
Return n indices in order
Cannot find a julia method which does, this.
If a faster method is found, replace it later.
"""
function indminn!(ivec, x::AbstractVector{Float64}, val::Float64)
	# using enumerate to avoid indexing
	n=length(ivec)
	if(length(x)<n)
		# if the desired indices is greater than length of x, what to do..
		fill!(ivec,one(eltype(ivec)))
	else
		fill!(ivec,zero(eltype(ivec)))
		for inn in 1:n
			ivec[inn] = indminimum(x, val, ivec)
		end
		sort!(ivec)
	end
	return ivec
end

"""
minimum index using imask
"""
function indminimum(x, val, imask=[])
	min_i=0
	min_x=typemax(Float64)
	for (i, xi) in enumerate(x)
		dist = abs(xi - val)
		if ((dist < min_x) && (i ∉ imask))
			min_x = dist
			min_i = i
		end
	end
	return min_i
end

function indminn(x::AbstractVector{Float64}, val::Float64, n::Int64=1)
	ivec=fill(0,n)
	indminn!(ivec,x,val)
	return ivec
end

"return index such that "
function indminn_inside(x, valbeg, valend)
	iibeg1 = indminimum(x, valbeg)
	iibeg2 = indminimum(x, valbeg,iibeg1)
	iiend1 = indminimum(x, valend)
	iiend2 = indminimum(x, valend,iiend1)

	if((valbeg-x[iibeg1])*(valend-x[iibeg1]) <= 0.0)
		iibeg=iibeg1
	else
		iibeg=iibeg2
	end
	if((valbeg-x[iiend1])*(valend-x[iiend1])<=0.0)
		iiend=iiend1
	else
		iiend=iiend2
	end

	return iibeg, iiend
end


