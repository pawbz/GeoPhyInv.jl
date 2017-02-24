module Data

using SeismicInversion.Acquisition
type DataTime
	δt :: Float64
	tvec :: Float64
	attrib::AbstractString
end


end # module
