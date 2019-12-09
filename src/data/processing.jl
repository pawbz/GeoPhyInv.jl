

"""
Normalize time-domain seismic data.

# Arguments 

* `data::TD` : input data
* `attrib::Symbol` : decide kind of normalization
  * `=:recrms` the record at every receiver is normalized with its RMS value
  * `=:recmax` the record at every receiver is normalized with its maximum value

# Return

* normalized data as `TD`
"""
function TD_normalize(data::TD, attrib::Symbol=:recrms)
	datan = deepcopy(data);
	TD_normalize!(datan, attrib)
	return datan
end
function TD_normalize!(data::TD, attrib::Symbol=:recrms)
	nr = data.ageom.nr;	nss = data.ageom.nss;	nt = length(data.tgrid);
	for ifield = 1:length(data.fields), iss = 1:nss
		dd=data.d[iss, ifield]
		scs=vecnorm(dd,2)

		for ir = 1:nr[iss]
			ddv=view(dd, :, ir)

			if(attrib == :recrms)
				sc=vecnorm(ddv,2)
				rmul!(ddv,inv(sc))
			elseif(attrib == :recmax)
				sc=vecnorm(ddv,Inf)
				rmul!(ddv,inv(sc))
			elseif(attrib == :srcrms)
				rmul!(ddv,inv(scs))
			else
				error("invalid attrib")
			end
		end
	end
	return data
end


function TD_filter!(data::TD; fmin=nothing, fmax=nothing)

	if((fmin===nothing) && (fmax===nothing))
		return nothing
	end
	fs = inv(step(data.tgrid));
	designmethod = Butterworth(4);
	filtsource = Bandpass(fmin, fmax; fs=fs);

	nr = data.ageom.nr;	nss = data.ageom.nss;	nt = length(data.tgrid);
	for ifield = 1:length(data.fields), iss = 1:nss
		dd=data.d[iss, ifield]
		scs=vecnorm(dd,2)

		for ir = 1:nr[iss]
			ddv=view(dd, :, ir)

			filt!(ddv, digitalfilter(filtsource, designmethod), ddv)
		end
	end
	return data
end

"""
Construct TD using data at all the unique receiver positions
for all supersources.

* `d::Array{Float64}` : the data matrix ordered in order such that time-domain modelling schemes are fast, i.e., [irec,ifield,it,nss]

"""
function TD_urpos(d::Array{Float64}, 
		   fields::Vector{Symbol}, 
		   tgrid::StepRangeLen, 
		   acq::AGeom,
		   nur::Int64,
		   urpos::Tuple{Array{Float64,1},Array{Float64,1}
		  }
		   )
	dout = [zeros(length(tgrid),acq.nr[iss]) for iss=1:acq.nss, ifield=1:length(fields)] 

	for ifield=1:length(fields), iss=1:acq.nss, ir=1:acq.nr[iss]
		# find index in urpos
		irr=find([[urpos[1][i]-acq.rz[iss][ir],
		       urpos[2][i]-acq.rx[iss][ir]] == [0., 0.,] for i in 1:nur])

		dout[iss, ifield][:,ir] = d[irr[1],ifield, :,iss] 
	end

	return TD(dout, fields, tgrid, acq)

end

"""
Apply source and receiver coupling functions to TD.
Currently, only source filters are applied.

# Arguments

* `s::TD` : input data
* `r::TD` : input data
* `w::Coupling.TD` : input source and receiver filters
* `attrib::Symbol` : attribute to 
  * `=:s` to apply `w` to `r` and modify `s`
  * `=:r` to apply adjoint of `w` to `s` and modify `r`
  * `=:w` modify `w` using `r` and `s`

TODO: need to work on parallelization and speed up here
"""
function TDcoup!(
               s::TD,
	       r::TD,
	       w::Coupling.TD,
	       attrib::Symbol
	       )
	nr = r.ageom.nr;	nss = r.ageom.nss;	nt = length(r.tgrid);
	fields = (w.fields == r.fields == s.fields) ? w.fields : error("different fields")
	sv=zeros(length(s.tgrid))
	rv=zeros(length(r.tgrid))
	wv=zeros(length(w.tgridssf))
	for ifield = 1:length(fields), iss = 1:nss, ir = 1:nr[iss]
		# receiver coupling
	#	conv!(s.d[iss, ifield][:, ir],r.d[iss, ifield][:, ir],
#		 w.rf[iss, ifield][:,ir],:s)
		# source coupling
		sv=s.d[iss, ifield][:, ir]
		rv=r.d[iss, ifield][:, ir]
		wv=w.ssf[iss, ifield]
#		conv!(sv, rv, wv, attrib)
		s.d[iss, ifield][:, ir]=copy(sv)
		r.d[iss, ifield][:, ir]=copy(rv)
		w.ssf[iss, ifield][:]=copy(wv)
	end
end # TDcoup


