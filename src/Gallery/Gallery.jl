module Gallery

using Pkg
import GeoPhyInv.IO
import GeoPhyInv: Medium, update!
import GeoPhyInv.FWI
import GeoPhyInv.Utils
using Statistics

global marmousi_folder=joinpath(@__DIR__, "marmousi2")



"""
Gallery of source signals `Src`.

# Arguments 
* `attrib::Symbol` : 
* `nss::Int64=1` : number of supersources

# Outputs
* `attrib=:acou_homo1` : 
"""
function Src(attrib::Symbol, nss::Int64=1)
	if(attrib == :acou_homo1)
		tgrid=range(0.0,stop=2.0,length=1000)
		wav = Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25, )
		return Src_fixed(nss, 1, [:P], wav, tgrid)
	elseif(attrib == :acou_homo2)
		tgrid=range(0.0,stop=2.0,length=250)
		wav = Utils.Wavelets.ricker(3.0, tgrid, tpeak=0.3, )
		return Src_fixed(nss, 1, [:P], wav, tgrid)
	elseif(attrib == :vecacou_homo1)
		tgrid=range(0.0,stop=2.0,length=1000)
		wav = Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25, )
		return Src_fixed(nss, 1, [:P, :Vx, :Vz], wav, tgrid)
	end
end


for file in ["fwi.jl", "fdtd.jl"]   
	fn=joinpath(@__DIR__, file)
	include(fn)
	#Revise.track(@__MODULE__,fn)
end

		 	
end # module
