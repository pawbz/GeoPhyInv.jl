mutable struct P_misfit
	x::TD # modelled data
	y::TD # observed data
	w::TD # weights
	dJx::TD # gradient w.r.t. x
	ynorm::Float64 # normalize functional with energy of y
end

function P_misfit(x, y; w=nothing)

	if(w===nothing) 
		w=Data.TD_ones(y.fields,y.tgrid,y.acqgeom)
	end
	!(isapprox(w,y)) && error("weights have to be similar to y")
	!(isapprox(x,y)) && error("cannot measure LS misfit b/w dissimilar data")

	!(isequal(x.acqgeom, y.acqgeom)) && error("observed and modelled data should have same acqgeom")
	!(isequal(x.fields, y.fields)) && error("observed and modelled data should have same fields")

	dJx=TD_zeros(x)

	ynorm = dot(y, y)
	(ynorm == 0.0) && error("y cannot be zero")

	pa=P_misfit(x,y,w,dJx,ynorm)

	return pa
end


function func_grad!(pa::P_misfit, grad=nothing)

	# some aliases

	y=pa.y
	x=pa.x
	tgrid = x.tgrid;
	acq = x.acqgeom;
	fields = x.fields;
	nss = length(acq);


	J = 0.0;
	for ifield=1:length(fields), iss=1:nss
		xr=x.d[iss,ifield]
		yr=y.d[iss,ifield]
		wr=pa.w.d[iss,ifield]

		if(grad===nothing)
			JJ=Misfits.error_squared_euclidean!(nothing,xr,yr,wr,norm_flag=false)
		elseif(grad==:dJx)
			dJxr=pa.dJx.d[iss,ifield]
			JJ=Misfits.error_squared_euclidean!(dJxr,xr,yr,wr,norm_flag=false)
		else
			error("invalid grad")
		end
		J += JJ 
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
mutable struct P_misfit_ssf
	x::TD # modelled data
	y::TD # observed data
	pa::P_misfit
	dJx::TD # gradient w.r.t. x
	xr::TD # modelled data after resampling
	xrc::TD # modelled data after resampling and convolution with source coupling
	dJxr::TD # gradient w.r.t. xr
	dJxrc::TD # gradient w.r.t. xr
	dJx::TD # gradient w.r.t. x
	dJssf::Matrix{Vector{Float64}} # gradient w.r.t source filters
	ynorm::Float64 # normalize functional with energy of y
	coup::Coupling.TD
	paconvssf::Vector{FBD.P_conv{Float64,2,2,1}} # convolutional model for source filters
	paconvrf::FBD.P_conv{Float64,1} # convolutional model for receiver filters
	pacse::Matrix{FBD.P_misfit_xcorr}
	func::Matrix{Function} # function to compute misfit for every ss and fields
	painterp::Interpolation.Param{Float64}
	interp_flag # decide if interpolation is necessary or not
end

function P_misfit_ssf(x, y; w=nothing, coup=nothing, func_attrib=:cls)

	if(coup===nothing)
		 coup=Coupling.TD_delta(y.tgrid,[0.1,0.1],[0.0, 0.0], [:P], y.acqgeom)
	end


	paconvssf=[FBD.P_conv(ssize=[coup.tgridssf.nx], 
			dsize=[length(y.tgrid),y.acqgeom.nr[iss]], 
			gsize=[length(y.tgrid),y.acqgeom.nr[iss]], 
		      slags=coup.ssflags, 
		      dlags=[length(y.tgrid)-1, 0], 
		      glags=[length(y.tgrid)-1, 0]) for iss in 1:y.acqgeom.nss]

	dJssf=deepcopy(coup.ssf)

	pacls=P_misfit(TD_zeros(y), y; w=w)

	!(isequal(x.acqgeom, y.acqgeom)) && error("observed and modelled data should have same acqgeom")
	!(isequal(x.fields, y.fields)) && error("observed and modelled data should have same fields")

	xr=TD_zeros(y)
	xrc=TD_zeros(y)
	dJxr=TD_zeros(y)
	dJxrc=TD_zeros(y)
	dJx=TD_zeros(x)

	ynorm = dot(y, y)
	(ynorm == 0.0) && error("y cannot be zero")

	dJxc=zeros(length(y.tgrid))

	if(func_attrib==:cls)
		pacse=[FBD.P_misfit_xcorr(1, 1,y=zeros(1,1)) for i in 1:2, j=1:2] # dummy
		func=[(dJx,x)->Misfits.error_squared_euclidean!(dJx,x,y.d[iss,ifield],w.d[iss,ifield]) for iss in 1:y.acqgeom.nss, ifield=1:length(y.fields)]
	elseif(func_attrib==:xcorrcls)
		pacse=[FBD.P_misfit_xcorr(length(y.tgrid), y.acqgeom.nr[iss],y=y.d[iss,ifield]) for iss in 1:y.acqgeom.nss, ifield=1:length(y.fields)]
		func=[(dJx,x)->FBD.func_grad!(dJx,x,pacse[iss,ifield]) for iss in 1:y.acqgeom.nss, ifield=1:length(y.fields)]
	end

	pa=P_misfit(x,y,w,xr,xrc,dJxr,dJxrc,
		    dJx,dJssf,ynorm,coup,paconvssf, paconvrf, pacse, func, painterp, interp_flag)


	return pa
end

function func_grad!(pa::P_misfit_ssf, grad=nothing)

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
