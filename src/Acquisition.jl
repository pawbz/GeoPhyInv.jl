module Acquisition

import SeismicInversion.Grid: Tim
import SeismicInversion.Wavelets

type Geom
	sx::Array{Float64}
	sz::Array{Float64}
	rx::Array{Float64}
	rz::Array{Float64}
	ns::Int64
	nr::Int64
end

function Geom()
	return Geom(1, 1)
end

function Geom(
	      ns::Int64,
	      nr::Int64
	     )
	sx = fill(-300.0,1)
	isequal(ns,1) ? sz = fill(-300.0,1) : sz = linspace(-300,300,ns)

	rx = fill(300.0,1)
	isequal(nr,1) ? rz = fill(-300.0,1) : rz = linspace(-300,300,nr)
	return Geom(sx, sz, rx, rz, ns, nr)
end


type Src
	wav::Array{Float64}
	tgrid::Tim
end

function Src()
	tgrid = Tim()
	wav = Wavelets.ricker()
	return Src(wav, tgrid)
end


end # module
