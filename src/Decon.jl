# blind deconvolution
__precompile__()

module Decon
import JuMIT.DSP
import JuMIT.Misfits
using Optim, LineSearches

type Param
	gf::Array{Float64,2}
	ntgf::Int64
	dobs::Array{Float64,2}
	dcal::Array{Float64,2}
	wav::Vector{Float64}
	nt::Int64
	nr::Int64
	attrib_inv::Symbol
	verbose::Bool
end

function Param(ntgf, nt, nr, gfobs, wavobs; verbose=false, attrib_inv=:gf,) 
	dobs=zeros(nt, nr)
	dcal=zeros(nt, nr)
	
	# initial values are random
	wav=randn(nt)
	gf=randn(ntgf,nr)

	pa = Param(gf, ntgf, dobs, dcal, wav, nt, nr, attrib_inv, verbose)

	forward_simulation!(pa, mattrib=:gfwav, dattrib=:dobs, gf=gfobs, wav=wavobs)

	return pa
end
function ninv(pa)
	if(pa.attrib_inv == :wav)
		return pa.nt
	else(pa.attrib_inv == :gf)
		return pa.ntgf*pa.nr
	end
end

function model_to_x!(x, pa)
	if(pa.attrib_inv == :wav)
		x[:] = copy(pa.wav)
	else(pa.attrib_inv == :gf)
		x[:] = vec(copy(pa.gf))
	end
	return x
end


function x_to_model!(x, pa)
	if(pa.attrib_inv == :wav)
		pa.wav[:] = copy(x)
	else(pa.attrib_inv == :gf)
		pa.gf = reshape(x, pa.ntgf, pa.nr)
	end
	return pa
end

function forward_simulation!(
			pa::Param,
			x::Vector{Float64}=zeros(ninv(pa)),
			last_x::Vector{Float64}=ones(ninv(pa));
			gf::Array{Float64}=zeros(pa.ntgf, pa.nr),
			wav::Array{Float64}=zeros(pa.nt),
			mattrib::Symbol=:x,
			dattrib::Symbol=:dcal,
			   )

	if(x!=last_x)

		if(mattrib == :x)
			x_to_model!(x, pa)
			gf=pa.gf
			wav=pa.wav
		elseif(mattrib == :gf)
			wav=pa.wav
			(iszero(gf)) && error("need gf")
		elseif(mattrib == :wav)
			gf=pa.gf
			(iszero(gf)) && error("need gf")
		elseif(mattrib == :gfwav)
			(iszero(gf)) && error("need gf")
			(iszero(wav)) && error("need wav")
		else
			error("invalid mattrib")
		end

		rv = zeros(pa.ntgf)
		sv = zeros(pa.nt)
		wv = zeros(pa.nt)
		#pa.verbose && println("updating buffer")
		dtemp = zeros(pa.dobs)
		last_x[:] = x[:]
		for ir in 1:pa.nr
			rv[:] = gf[:, ir]
			wv[:] = wav
			DSP.fast_filt!(sv, rv, wv, :s, nwplags=pa.nt-1)
			dtemp[:, ir] = copy(sv)
		end

		setfield!(pa, dattrib, deepcopy(dtemp))
		return pa
	end
end

function func_grad!(storage, x::Vector{Float64},
			last_x::Vector{Float64}, 
			pa)

	#pa.verbose && println("computing gradient...")

	# x to w 
	x_to_model!(x, pa)

	rv = zeros(pa.ntgf)
	sv = zeros(pa.nt)
	wv = zeros(pa.nt)
	forward_simulation!(pa, x, last_x)

	if(storage === nothing)
		# compute misfit and δdcal
		f = Misfits.fg_cls!(nothing, pa.dcal, pa.dobs)
		return f
	else
		ddcal = similar(pa.dobs)
		f = Misfits.fg_cls!(ddcal, pa.dcal, pa.dobs)

		gwav = similar(pa.wav)
		gwav[:] = 0.0
		ggf = similar(pa.gf)
		if(pa.attrib_inv == :wav)
			for ir in 1:pa.nr
				rv[:] = pa.gf[:, ir]
				sv[:] = ddcal[:, ir]
				DSP.fast_filt!(sv, rv, wv, :w, nwplags=pa.nt-1)
				gwav[:] += wv
			end
			storage[:] = vec(gwav)
		else(pa.attrib_inv == :gf)
			for ir in 1:pa.nr
				sv = ddcal[:, ir]
				wv = pa.wav
				DSP.fast_filt!(sv, rv, wv, :r, nwplags=pa.nt-1)
				ggf[:, ir] = copy(rv)
			end
			storage[:] = copy(vec(ggf))
			return storage
		end

	end
end



# core algorithm
function update!(pa::Param; store_trace::Bool=true, extended_trace::Bool=false, time_limit=Float64=2.0*60., 
	       f_tol::Float64=1e-8, g_tol::Float64=1e-8, x_tol::Float64=1e-8)

	# convert initial model to the inversion variable
	x = zeros(ninv(pa));
	last_x = rand(size(x)) # reset last_x

	#pa.verbose && println("updating pa...")
	#pa.verbose && println("> number of inversion variables:\t", length(x)) 

	# initial w to x
	model_to_x!(x, pa)

	df = OnceDifferentiable(x -> func_grad!(nothing, x, last_x, pa),
		  (storage, x) -> func_grad!(storage, x, last_x, pa))
	"""
	Unbounded LBFGS inversion, only for testing
	"""
	res = optimize(df, x, 
		       LBFGS(),
		       Optim.Options(g_tol = g_tol, f_tol=f_tol, x_tol=x_tol,
		       iterations = 100, store_trace = store_trace,
		       extended_trace=extended_trace, show_trace = false))
	pa.verbose && println(res)

	x_to_model!(Optim.minimizer(res), pa)

	return res

end


function update_all!(pa)

	converged_all=false
	itr_all=0
	maxitr_all=10
	while !converged_all && itr_all < maxitr_all
		itr_all += 1
		pa.verbose && (itr_all > 1) && println("\t", "failed to converge.. reintializing round trips")
		# starting models
		pa.wav[:] = randn(size(pa.wav))
		pa.gf[:] = randn(size(pa.gf))
	
		fgfmin=Inf
		fgfmax=0.0
		fwavmin=Inf
		fwavmax=0.0

		maxitr=200
		itr=0
		converged=false
		tol=1e-8
		while !converged && itr < maxitr
			pa.verbose && (itr > 1) && println("\t", "round trip: ", itr)
			itr += 1
			pa.attrib_inv=:gf    
			resgf = update!(pa)
			fgf = Optim.minimum(resgf)
			fgfmin = min(fgfmin, fgf) 
			fgfmax = max(fgfmax, fgf) 
			pa.attrib_inv=:wav    
			reswav = update!(pa)
			fwav = Optim.minimum(reswav)
			fwavmin = min(fwavmin, fwav) 
			fwavmax = max(fwavmax, fwav) 

			converged = ((abs(fwav-fgf)) < tol)
		end
		converged_all = (fgfmin/fgfmax < 1e-3) & (fwavmin/fwavmax < 1e-3)
	end
end




# initialize gf and ssf 

# create unknown vector

end # module
