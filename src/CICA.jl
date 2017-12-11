__precompile__()


"""
Independent Component Analysis for Convolutive Mixtures
"""
module CICA

#using MultivariateStats
using Distributions
import JuMIT.Grid:M1D
import JuMIT.Misfits
import JuMIT.Data
import JuMIT.Grid
import JuMIT.DSP
import JuMIT.Acquisition

type ICA{T}
	"blended data at receivers(nr, nt)"
	x::Array{T,2}
	"blended data after whitening (ns, nt)"
	xhat::Array{T,2}
	"deblended data (ns, nt)"
	s::Array{T,2}
	"number of sources or independent components"
	ns::Int64
	"number of receivers"
	nr::Int64
	"number of samples"
	nt::Int64
	"store G(Wx)"
	GWx::Matrix{T}
	"store g(Wx)"
	gWx::Matrix{T}
	"store dg(Wx)"
	dgWx::Matrix{T}
	"whitening matrix --- this is updated after preprocessing"
	Q::Array{T,3}
	"unmixing matrix --- this is updated after ICA"
	W::Array{T,3}
	"previous unmixing matrix"
	Wp::Array{T,3}
	"store mean vector"
	mv::Array{T,2}
	"real or complex ica"
	attrib::Symbol
	"maximum iterations"
	maxiter::Int64
	"receiver at which deblending happens"
	magic_recv::Int64
	"tolerance"
	tol::Float64
	"bins"
	bins::Vector{UnitRange{Int64}}
end


function ICA(dat::Data.TD; 
	     ica_func::Function=(x,ns)->ICA(x,ns),)
	Acquisition.Geom_isfixed(dat.acqgeom) || error("implemented only for fixed spread geom")
	acqgeom=dat.acqgeom
	ns=acqgeom.nss
	nr=acqgeom.nr[1] # fixed spread
	nt=dat.tgrid.nx
	xx=zeros(nt, nr, ns)
	for iss in 1:ns
		for ir in 1:nr
			for it in 1:nt
				xx[it,ir,iss] = dat.d[iss,1][it,ir]
			end
		end
	end
	# real fft gves only +ve freqs
	xx=rfft(xx, [1])

	xstack=complex.(zeros(size(xx,2),size(xx,1)))
	A=zeros(nr, ns)
	for iss in 1:ns
		for ir in 1:nr
			for i in 1:size(xx,1)
				# stack over supersources
				xstack[ir,i] += xx[i,ir,iss]
			end
			# mixing matrix (wrong!!!)
			A[ir,iss]=var(xx[:,ir,iss])
		end
	end

	ica = ica_func(xstack, ns), complex.(A)
	ica=ica[1]

	fastica!(ica)

	s=transpose(ica.s)

	xx=irfft(s, nt, [1])


	return xx

end
"""
`ns`	: 
`dat`	: dat


"""
function ICA(x,
	     ns;
	     maxiter::Int=1000,
	     tol::Real=1e-10,                 # convergence tolerance
	     Win=nothing,  # init guess of W, size (m, k)
	     nbins=1, # number of bins
	     )

	T = eltype(x)
	nr, nt = size(x)

	# check input arguments
	(nt < 1) && error("there must be more than one samples, i.e. nt > 1.")
	(ns > nr) && error("ns must not exceed nr")

	(maxiter < 1) && error("maxiter must be greater than 1.")
	(tol < 0.) && error("tol must be positive.")


	# bins
	bindices=[round(Int, s) for s in linspace(0,nt,nbins+1)]
	bins=Array{UnitRange{Int64}}(nbins)
	for ib in 1:nbins           
		bins[ib]=bindices[ib]+1:bindices[ib+1]
	end

	println("FastICA\t")
	println("===================")
	println("total samples\t",nt)
	println("total number of bins\t",nbins)
	println("average samples in each bin\t",round(Int,mean(length.(bins))))
	println("")

	# allocating vectors, that store G, g, and dg
	GWx = Matrix{T}(ns, nt)
	gWx = Matrix{T}(ns, nt)
	dgWx = Matrix{T}(ns, nt)

	# random initialization of unmixing matrix
	W = Array{T}(ns,ns, nbins)
	Wp = Array{T}(ns,ns, nbins)


	Q = Array{T}(nr, ns, nbins)
	xhat = Matrix{T}(ns, nt)

	mv = Matrix{T}(nr, nbins)
	# deblended data
	s = Matrix{T}(ns, nt)

	(T == Float64) ? attrib=:real : attrib=:complex

	ica=ICA(x, xhat, s, ns, nr, nt, GWx, gWx, dgWx, Q, W, Wp, mv, attrib,
	 maxiter, 1, tol, bins)

	unmixing_initialize!(ica, Win)

	return ica
end

function unmixing_initialize!(ica::ICA, Win=nothing)
	T = eltype(ica.x)
	ns=ica.ns
	W=ica.W
	Wp=ica.Wp
	bins=ica.bins

	if(!(Win===nothing))
		copy!(W, Win)
	else
		if(ica.attrib==:real)
			for ib in 1:length(bins), is1=1:ns, is2=1:ns,
				W[is1,is2,ib]=randn()
			end
		elseif(ica.attrib==:complex)
			for ib in 1:length(bins), is1=1:ns, is2=1:ns,
				W[is1,is2,ib]=complex.(randn()+randn())
			end
		else
			error("unknown attrib")
		end
	end

	# normalize each column of W necessary if arbitary Win is given
	for ib in 1:length(bins)
		for is = 1:ns
			w = view(W,:,is,ib)
			scale!(w, 1.0 / sqrt(sum(abs2, w)))
		end
	end
	return ica
end


"""
* update mv by storing mean of each component of x
* remove mean from each component of x
"""
function remove_mean!(ica::ICA)
	for ib in 1:length(ica.bins)
		for ir=1:ica.nr
			ica.mv[ir, ib]=zero(eltype(ica.mv))
			# compute mean
			for it=ica.bins[ib]
				ica.mv[ir,ib] += ica.x[ir,it]
			end
			ica.mv[ir, ib] *= inv(length(ica.bins[ib]))
			# remove mean
			for it=ica.bins[ib]
				ica.x[ir,it] -= ica.mv[ir,ib]

			end
		end
	end
end

"""
add stored mean to x
"""
function add_mean!(ica::ICA)
	for ib in 1:length(ica.bins)
		for ir=1:ica.nr
			# add mean
			for it=ica.bins[ib]
				ica.x[ir, it] += ica.mv[ir, ib]
			end
			ica.mv[ir, ib] = zero(eltype(ica.mv))
		end
	end
end

"""
Apply deblending matrix W to preprocessed data xhat and return deblended data ica.s
"""
function deblend!(ica::ICA)
	for ib in 1:length(ica.bins)
		s=view(ica.s,:,ica.bins[ib])
		W=view(ica.W,:,:,ib)
		xhat=view(ica.xhat,:,ica.bins[ib])
		Ac_mul_B!(s, W, xhat)
	end
	return ica.s
end

"""
update xhat
"""
function preprocess!(ica::ICA)
	# remove and store mean
	remove_mean!(ica)

	do_whiten!(ica)
	 
	# adding mean back to x; avoided making a copy of x
	add_mean!(ica)
end

# Fast ICA
#
# Reference:
#
#   Aapo Hyvarinen and Erkki Oja
#   Independent Component Analysis: Algorithms and Applications.
#   Neural Network 13(4-5), 2000.
#
function fastica!(ica::ICA; verbose::Bool=false, A=nothing)   

	preprocess!(ica)


	ns=ica.ns
	nt=ica.nt
	bins=ica.bins
	nbins=length(ica.bins)

	# aliases for x, W, Wp
	x=ica.xhat
	s=ica.s
	W=ica.W
	Q=ica.Q
	Wp=ica.Wp

	# some other aliases
	maxiter = ica.maxiter

	# aliases GWx, gWx, dgWx
	GWx = ica.GWx;	gWx = ica.gWx;	dgWx = ica.dgWx

	# main loop
	t=0
	converged=false
	while !converged && t < maxiter
		t += 1
		copy!(Wp, W)
		deblend!(ica)

		ICAfunc!(s, GWx, gWx, dgWx, ica.attrib)

		updateW!(W,gWx,dgWx,s,x,bins,ica.attrib)

		symmetric_decorrelation!(W)

		converged=compareWs(W,Wp,ica.tol)


		# compute the distance between A and W'Q
		# this gives the measure of the success of the algorithm
		# note that W'Q and A can still be different by a 
		# permutation and scaling matrix, so...
		# if(!(A === nothing))
		#	SSE[t] = error_scale_perm(W'*Q', A)
		# end
	end

	# store expectation of G -- it measures Gaussianity
	EG = zeros(ns, nbins)

	# storing the expectation of G(Wx) --- measure of Gaussianity
	for ib=1:nbins
		for it in bins[ib]
			for ic=1:ns
				EG[ic,ib] += (GWx[ic, it])
			end
		end
		for ic=1:ns
			EG[ic,ib] /= length(bins[ib])
		end
	end

	fix_permutation!(W, EG)

	fix_scaling!(ica)

	# store error in the unmixing matrix -- only for testing (commented)
	SSE = (!(A === nothing)) ? error_unmixing(ica, A) : zeros(nbins)

	# perform final deblending after fixing permutation and scaling
	deblend!(ica)

	# test the performance of the algorithm
	return ica, SSE, EG
end

function error_unmixing(ica, A)        
	nbins=length(ica.bins)     
	Q=ica.Q        
	W=ica.W        
	SSE=zeros(nbins)
	for ib=1:nbins                 
		QQ=view(Q,:,:,ib)                 
		WW=view(W,:,:,ib)                 
		SSE[ib] = error_scale_perm(WW'*QQ', A)        
	end
	return SSE
end

"""
compute the distance between s and 
* `s` actual source matrix [nt, ns]
* `ica` ica output s is used

There is still global permutation 
and scaling to be fixed.
"""
function error_after_ica(ica, s)
	(size(s) ≠ size(ica.s)) && error("dimension")
	
	errr=zeros(size(s,1))

	for is in 1:2
		ss0=view(s,is,:)
		err=zeros(size(s,1))
		for is1 in 1:2
			ss=view(ica.s,is1,:)
			err[is1] = Misfits.error_after_scaling(ss,ss0)[1]
		end
		errr[is] = minimum(err)
	end
	return errr
end

"""
Minimum Distortion Principle.
This method is not yet memory optimized.
Performs ICA for convolutive mixtures.
# Arguments
# Output
* deblended data at `magic_recv`
"""
function fix_scaling!(ica)
	nbins=length(ica.bins)
	ns=ica.ns
	W=ica.W
	Q=ica.Q
	WQ=zeros(eltype(W),size(W,1),size(W,2))
	for ib in 1:nbins
		WW=view(W,:,:,ib)
		QQ=view(Q,:,:,ib)

		Ac_mul_Bc!(WQ, WW, QQ)
		# A is the mixing matrix
		A = inv(WQ);

		# normalize mixing matrix, w.r.t. magic_recv
		A[:,1] /= A[ica.magic_recv,1];
		A[:,2] /= A[ica.magic_recv,2];

	 	copy!(WQ, (inv(A)));
		
		# scalings are only put in WW, so using inv(QQ)
		copy!(WW, inv(QQ)*WQ')
	end
end

"""
Order the components according to the EG values.
EG measures Gaussianity and hence, we are ordering sources 
in each bin according to their distance from Gaussian.
"""
function fix_permutation!(W, EG)
	nbins=size(W,3)
	for ib=1:nbins
		WW=view(W,:,:,ib)
		# permutation order
		ord = sortperm(EG[:,ib]; rev=true)
		copy!(WW, WW[:,ord])
	end
end

function updateW!(W,gWx,dgWx,s,x,bins,attrib)
	nbins=length(bins)
	ns=size(W,2)
	if(attrib == :complex)

		# loop over bins; W for each bin is updated
		for ib=1:nbins
			for ic=1:ns
				b=zero(eltype(W))
				@simd for it in bins[ib]
					b+=gWx[ic,it]+(abs.(s[ic,it]).^2*dgWx[ic,it])
				end
				b /= length(bins[ib])
				
				for ic1=1:ns
					a=zero(eltype(x))
					@simd for it in bins[ib]
						a+=x[ic1,it]*conj(s[ic,it])*gWx[ic,it]
					end
					a /= length(bins[ib])
					W[ic1, ic, ib]=a-b*W[ic1,ic, ib]
				end
			end
		end
	elseif(attrib == :real)
		# loop over bins; W for each bin is updated
		for ib=1:nbins
			for ic = 1:ns
				b=zero(eltype(x))
				@simd for it=bins[ib]
					b+=dgWx[ic,it]
				end
				b /= length(bins[ib])
				for ic1=1:ns
					a=zero(eltype(x))
					@simd for it=bins[ib]
						a+=x[ic1,it]*gWx[ic,it]
					end
					a /= length(bins[ib])
					W[ic1, ic, ib]=a-b*W[ic1,ic, ib]
				end
			end
		end
	else
		error("invalid attrib")
	end
end

function ICAfunc!(s, GWx, gWx, dgWx, attrib)

	if(attrib==:complex)
		# for complex ICA (they used this constant in the paper)
		eps = 0.1
		@simd for it=1:size(s,2)
			for is=1:size(s,1)
				ss=eps+abs(s[is, it])^2
				GWx[is, it] = log(ss)
				gWx[is, it] = inv(ss)
				dgWx[is, it] = -1. * inv(ss*ss)
			end
		end
	elseif(attrib==:real)
		@simd for it=1:size(s,2)
			for is=1:size(s,1)
				GWx[is,it]=exp(-0.5*s[is,it]*s[is,it])
				gWx[is,it]=s[is,it]*GWx[is,it]
				dgWx[is,it]=GWx[is,it]*(1.-(s[is,it]*s[is,it]))
			end
		end
	end
end


# compare W with Wp as a convergence criteria
function compareWs(W,Wp,tol)
	nbins=size(W,3)
	ns=size(W,2)
	chgall=zeros(nbins)
	for ib in 1:nbins
		for ic in 1:ns
			chg=0.0
			for ic1 in 1:ns
				chg += abs(W[ic1, ic, ib]-Wp[ic1, ic, ib]) 
			end
			chgall[ib] = max(chg, chgall[ib]) 
		end
	end
	return all(chgall .< tol)
end

function symmetric_decorrelation!(W)
	# Symmetric decorrelation
	for ib=1:size(W,3)
		WW=view(W,:,:,ib)
		copy!(WW, WW * sqrtm(inv(WW'*WW)))
	end
end

"""
Compute the error in the unmixing matrix.

* x is the unmixing matrix
* y is the original mixing matrix 

Ideally, x should be inverse of y upto scaling and permutation.
}	
"""
function error_scale_perm{T}(x::Matrix{T}, y::Matrix{T}) 
	K=abs.(x*y);
	nc=size(x,2)
	err=0.0
	for ic in 1:nc    
		kk = K[ic,:]    
		ii = indmax(abs.(kk))    
		kk ./= kk[ii]    
		err += sum(abs.(kk[1:ii-1])) + sum(abs.(kk[ii+1:end]))
	end
	return err
end



"""
Whiten the data that are input to ICA. Whietning matrix is stored in `Q`.
"""
function do_whiten!(ica::ICA)
    nr = ica.nr
    for ib in 1:length(ica.bins)
	    x=view(ica.x,:,ica.bins[ib])
	    xhat=view(ica.xhat,:,ica.bins[ib])
	    Q=view(ica.Q,:,:,ib)
	    C=Matrix{eltype(x)}(nr, nr)
	    ntI = 1/length(ica.bins[ib])

	    # construct covariance matrix
	    A_mul_Bc!(C, x, x)
	    scale!(C, ntI)

	    # eigenvalue decomposition of the covariance matrix
	    Efac = eigfact(C)
	    ord = sortperm(Efac.values; rev=true)
	    v, P = Efac.values[ord], Efac.vectors[:, ord]
	    # whitening matrix
	    copy!(Q, scale!(P, sqrt.(inv.((v)))))

	    # for testing --  this should be a diagonal matrix
	    #println(Q' * C * Q)
	    Ac_mul_B!(xhat, Q, x)
    end
    return ica.Q
end



end # module
