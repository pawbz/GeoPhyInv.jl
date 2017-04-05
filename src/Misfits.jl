module Misfits



"""
Normalized least-squares error between two vectors after 
estimating a scalar that best fits on to another.
Return misfit and α such that αx-y is minimum.
Normalization is done with respect to the 
norm of y.
"""
function error_after_scaling(
			     x::Array{Float64},
			     y::Array{Float64}
			    )

length(x) != length(y) ? error("vectors of different length") :

α = dot(x,y)/dot(x,x)

J = norm(y-α*x)/norm(y)
#J = norm(α.*x-y)/norm(y)

return J, α

end



end # module
