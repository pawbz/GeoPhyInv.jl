

"""
Apply different weighting functions to `TD`. Use this method to create a 
data preconditioning matrix

# Arguments Modified

* `dw::TD` : 

# Keyword Arguments

* `offsetlim::Vector{Float64}=[-Inf,Inf]` : [xoffsetlim, zoffsetlim], where the records with offsets > `offsetlim` are given zero weight 
* `tlim::Vector{Float64}=[dw.tgrid.x[1], dw.tgrid.x[end]]` : [tminimum, tmaximum], time mute window 
* `offsetpow::Vector{Float64}=[0.0,0.0]` : 
* `tpow::Float64=0.0` :
* `ttaperperc::Float64=0.` : taper window percentage for time window

* NOTE: if more than one simultaneous source are present, their mean position is considered to calculate offset.
"""
function TD_weight!(
		    dw::TD;
		    offsetlim::Vector{Float64}=[-Inf,Inf],
		    tlim::Vector{Float64}=[dw.tgrid.x[1], dw.tgrid.x[end]],
		    offsetpow::Vector{Float64}=[0.0,0.0],
		    tpow::Float64=0.0,
		    ttaperperc::Float64=0.,
		  )

	tvecexp = dw.tgrid.x - dw.tgrid.x[1]
	tmaxI = maximum(inv.(abs.(tvecexp)))
	nt = dw.tgrid.nx
	fields=dw.fields
	nss=dw.acqgeom.nss
	rx=dw.acqgeom.rx; rz=dw.acqgeom.rz
	sx=dw.acqgeom.sx; sz=dw.acqgeom.sz
	nr=dw.acqgeom.nr; ns=dw.acqgeom.ns

	itlim = sort(broadcast(indmin,[abs.(dw.tgrid.x-tlim[i]) for i in 1:2]))
	twin=zeros(nt)
	twin[itlim[1] : itlim[2]] = Signals.DSP.taper(ones(itlim[2]-itlim[1]+1),ttaperperc) 

	for ifield = 1:length(fields), iss = 1:nss
		zo = sqrt.((rz[iss][:]-mean(sz[iss])).^2) # offsets computed using mean of source position
		xo = sqrt.((rx[iss][:]-mean(sx[iss])).^2)
		zomaxI = maximum(zo)^(-1)
		xomaxI = maximum(xo)^(-1)
		for ir = 1:nr[iss]
			inoffsetlim = ( (abs.(zo[ir]) < abs.(offsetlim[2])) & (abs.(xo[ir]) < abs.(offsetlim[1])) )
			#if(inoffsetlim)
			#	# apply xoffset weighting 
			#	if(!isinf(xomaxI))
			#		dw.d[iss,ifield][:,ir] .*= exp(xo[ir]*offsetpow[1]*xomaxI) 
			#	end
			#	# apply zoffset weighting
			#	if(!isinf(zomaxI))
			#		dw.d[iss,ifield][:,ir] .*= exp(zo[ir]*offsetpow[2]*zomaxI)
			#	end
			#else
			#	dw.d[iss,ifield][:,ir] = 0.
			#end
			for it=1:nt
				dw.d[iss,ifield][it,ir] *= twin[it] #* exp(tvecexp[it]*tpow*tmaxI)
			end
		end
	end

end


