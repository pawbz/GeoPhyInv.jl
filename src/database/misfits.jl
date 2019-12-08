mutable struct VNamedD_misfit{T}
	x::Vector{NamedD{T}} # modelled data
	y::Vector{NamedD{T}} # modelled data
	w::Vector{NamedD{T}} # modelled data
	dJx::Vector{NamedD{T}} # modelled data
	ynorm::Float64 # normalize functional with energy of y
end

function VNamedD_misfit(x, y; w=nothing)

	if(w===nothing) 
		w=deepcopy(y);
		fill!(w, 1.0)
	end
	!(issimilar(w,y)) && error("weights have to be similar to y")
	!(issimilar(x,y)) && error("cannot measure LS misfit b/w dissimilar data")

	dJx=deepcopy(x); fill!(dJx,0.0)

	ynorm = dot(y, y)
	(ynorm == 0.0) && error("y cannot be zero")

	pa=VNamedD_misfit(x,y,w,dJx,ynorm)

	return pa
end


function func_grad!(pa::VNamedD_misfit, grad=nothing)

	# some aliases

	y=pa.y
	x=pa.x

	J = 0.0;
	for iss in 1:length(x)
		for ifield=names(x[iss].d)[1]
			xr=x[iss].d[ifield]
			yr=y[iss].d[ifield]
			wr=pa.w[iss].d[ifield]

			if(grad===nothing)
				JJ=Misfits.error_squared_euclidean!(nothing,xr,yr,wr,norm_flag=false)
			elseif(grad==:dJx)
				dJxr=pa.dJx[iss].d[ifield]
				JJ=Misfits.error_squared_euclidean!(dJxr,xr,yr,wr,norm_flag=false)
			else
				error("invalid grad")
			end
			J += JJ 
		end
	end
	return J
end


#=

"""
Calculate the distance between the observed data `y` and the calculated data `x`.
The time grid of the observed data can be different from that of the modelled data.
The acqistion geometry of both the data sets should be the same.

If `J` is the distance, the gradient of the misfit w.r.t to the calculated data is returned as `dJx`
* `w` used for data preconditioning
* `coupling` source and receiver coupling functions
"""
mutable struct VNamedD_misfit_ssf
	x::VNamedD # modelled data
	y::VNamedD # observed data
	pa::VNamedD_misfit
	dJx::VNamedD # gradient w.r.t. x
	xr::VNamedD # modelled data after resampling
	xrc::VNamedD # modelled data after resampling and convolution with source coupling
	dJxr::VNamedD # gradient w.r.t. xr
	dJxrc::VNamedD # gradient w.r.t. xr
	dJx::VNamedD # gradient w.r.t. x
	dJssf::Matrix{Vector{Float64}} # gradient w.r.t source filters
	ynorm::Float64 # normalize functional with energy of y
	coup::Coupling.VNamedD
	paconvssf::Vector{FBD.P_conv{Float64,2,2,1}} # convolutional model for source filters
	paconvrf::FBD.P_conv{Float64,1} # convolutional model for receiver filters
	pacse::Matrix{FBD.VNamedD_misfit_xcorr}
	func::Matrix{Function} # function to compute misfit for every ss and fields
	painterp::Interpolation.Param{Float64}
	interp_flag # decide if interpolation is necessary or not
end

function VNamedD_misfit_ssf(x, y; w=nothing, coup=nothing, func_attrib=:cls)

	if(coup===nothing)
		 coup=Coupling.VNamedD_delta(y.tgrid,[0.1,0.1],[0.0, 0.0], [:P], y.acqgeom)
	end


	paconvssf=[FBD.P_conv(ssize=[coup.tgridssf.nx], 
			dsize=[length(y.tgrid),y.acqgeom.nr[iss]], 
			gsize=[length(y.tgrid),y.acqgeom.nr[iss]], 
		      slags=coup.ssflags, 
		      dlags=[length(y.tgrid)-1, 0], 
		      glags=[length(y.tgrid)-1, 0]) for iss in 1:y.acqgeom.nss]

	dJssf=deepcopy(coup.ssf)

	pacls=VNamedD_misfit(VNamedD_zeros(y), y; w=w)

	!(isequal(x.acqgeom, y.acqgeom)) && error("observed and modelled data should have same acqgeom")
	!(isequal(x.fields, y.fields)) && error("observed and modelled data should have same fields")

	xr=VNamedD_zeros(y)
	xrc=VNamedD_zeros(y)
	dJxr=VNamedD_zeros(y)
	dJxrc=VNamedD_zeros(y)
	dJx=VNamedD_zeros(x)

	ynorm = dot(y, y)
	(ynorm == 0.0) && error("y cannot be zero")

	dJxc=zeros(length(y.tgrid))

	if(func_attrib==:cls)
		pacse=[FBD.VNamedD_misfit_xcorr(1, 1,y=zeros(1,1)) for i in 1:2, j=1:2] # dummy
		func=[(dJx,x)->Misfits.error_squared_euclidean!(dJx,x,y.d[iss,ifield],w.d[iss,ifield]) for iss in 1:y.acqgeom.nss, ifield=1:length(y.fields)]
	elseif(func_attrib==:xcorrcls)
		pacse=[FBD.VNamedD_misfit_xcorr(length(y.tgrid), y.acqgeom.nr[iss],y=y.d[iss,ifield]) for iss in 1:y.acqgeom.nss, ifield=1:length(y.fields)]
		func=[(dJx,x)->FBD.func_grad!(dJx,x,pacse[iss,ifield]) for iss in 1:y.acqgeom.nss, ifield=1:length(y.fields)]
	end

	pa=VNamedD_misfit(x,y,w,xr,xrc,dJxr,dJxrc,
		    dJx,dJssf,ynorm,coup,paconvssf, paconvrf, pacse, func, painterp, interp_flag)


	return pa
end

function func_grad!(pa::VNamedD_misfit_ssf, grad=nothing)

	tgrid = pa.x.tgrid;
	acq = pa.x.acqgeom;
	fields = pa.x.fields;
	nss = acq.nss;

	if(pa.interp_flag)
		# resample xr <-- x
		interp_spray!(pa.x, pa.xr, :interp, pa=pa.painterp)
		println("ffffff, HELEL")
	else
		copyto!(pa.xr,pa.x)
	end

	xr=pa.xr
	xrc=pa.xrc
	y=pa.y


	J = 0.0;
	for ifield=1:length(fields), iss=1:acq.nss
		xrr=xr.d[iss,ifield]
		xrcc=xrc.d[iss,ifield]
		wav=pa.coup.ssf[iss,ifield]

		nt=size(xrr,1)
		nr=size(xrr,2)

		paconv=pa.paconvssf[iss]
		# xrc <- xr apply source filter to xr
		FBD.mod!(paconv, :d, g=xrr, s=wav, d=xrcc)

		if(grad===nothing)
			JJ=pa.func[iss,ifield](nothing,  xrcc)
		elseif((grad==:dJx) || (grad==:dJssf))
			dJxrcc=pa.dJxrc.d[iss,ifield]
			JJ=pa.func[iss,ifield](dJxrcc,  xrcc)
		else
			error("invalid grad")
		end

		# dJxr <- dJxrc  apply adjoint of source filter to dJxc
		if(grad==:dJx)
			dJxrr=pa.dJxr.d[iss,ifield]
			dJxrcc=pa.dJxrc.d[iss,ifield]
			FBD.mod!(paconv, :g, g=dJxrr, s=wav, d=dJxrcc)
		end

		# dJssf <- dJxrc 
		if(grad==:dJssf)
			dJwav=pa.dJssf[iss,ifield]
			FBD.mod!(paconv, :s, g=xrr, s=dJwav, d=dJxrcc)
		end
		J += JJ 
	end

	if(pa.interp_flag)
		# spray dJx <-- dJxr
		interp_spray!(pa.dJx, pa.dJxr, :spray, pa=pa.painterp)
		println("ffffff, HELEL")
	else
		copyto!(pa.dJx,pa.dJxr)
	end

	(J == 0.0) && warn("misfit computed is zero")

	return J

end

=#
