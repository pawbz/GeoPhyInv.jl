__precompile__()


"""
Independent Component Analysis for Convolutive Mixtures
"""
module CICA

#using MultivariateStats
using Distributions
import JuMIT.Grid:M1D
import JuMIT.Misfits

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
	Q::Array{T}
	"unmixing matrix --- this is updated after ICA"
	W::Array{T}
	"previous unmixing matrix"
	Wp::Array{T}
	"store mean vector"
	mv::Matrix{T}
	"real or complex ica"
	attrib::Symbol
	"maximum iterations"
	maxiter::Int64
	"tolerance"
	tol::Float64
end

function ICA(x,
	     ns,
	     maxiter::Int=1000,
	     tol::Real=1e-10,                 # convergence tolerance
	     Win=nothing,  # init guess of W, size (m, k)
	     )

	T = eltype(x)

	# check input arguments
	nr, nt = size(x)
	(nt < 1) && error("there must be more than one samples, i.e. nt > 1.")
	(ns > nr) && error("ns must not exceed nr")

	(maxiter < 1) && error("maxiter must be greater than 1.")
	(tol < 0.) && error("tol must be positive.")

	# allocating vectors, that store G, g, and dg
	GWx = Matrix{T}(ns, nt)
	gWx = Matrix{T}(ns, nt)
	dgWx = Matrix{T}(ns, nt)

	# random initialization of unmixing matrix
	W = Matrix{T}(ns,ns)
	Wp = Matrix{T}(ns,ns)


	Q = Matrix{T}(nr, ns)
	xhat = Matrix{T}(ns, nt)

	mv = Matrix{T}(nr, 1)
	# deblended data
	s = Matrix{T}(ns, nt)

	(T == Float64) ? attrib=:real : attrib=:complex

	ica=ICA(x, xhat, s, ns, nr, nt, GWx, gWx, dgWx, Q, W, Wp, mv, attrib, maxiter, tol)

	unmixing_initialize!(ica, Win)

	return ica
end

function unmixing_initialize!(ica::ICA, Win=nothing)
	T = eltype(ica.x)
	ns=ica.ns
	W=ica.W
	Wp=ica.Wp

	if(!(Win===nothing))
		copy!(W, Win)
	else
		if(ica.attrib==:real)
			copy!(W,randn(ns,ns)) 
		elseif(ica.attrib==:complex)
			copy!(W,complex.(randn(ns,ns)+randn(ns,ns)))
		else
			error("unknown attrib")
		end
	end

	# normalize each column of W necessary if arbitary Win is given
	for is = 1:ns
		w = view(W,:,is)
		scale!(w, 1.0 / sqrt(sum(abs2, w)))
	end
	return ica
end


"""
* update mv by storing mean of each component of x
* remove mean from each component of x
"""
function remove_mean!(ica::ICA)
	mean!(ica.mv, ica.x)
	for it=1:ica.nt, ir=1:ica.nr
		ica.x[ir, it] -= ica.mv[ir, 1]
	end
end

"""
add stored mean to x
"""
function add_mean!(ica::ICA)
	for it=1:ica.nt, ir=1:ica.nr
		ica.x[ir, it] += ica.mv[ir, 1]
	end
	for ir=1:ica.nr
		ica.mv[ir, 1] = zero(eltype(ica.mv))
	end
end

"""
Apply deblending matrix W to preprocessed data xhat and return deblended data ica.s
"""
function deblend!(ica::ICA)
	Ac_mul_B!(ica.s, ica.W, ica.xhat)
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

	ns, nt = ica.ns, ica.nt
	ntI = inv(nt)

	preprocess!(ica)

	# for complex ICA
	eps = 0.1

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

	# store expectation of G -- it measures Gaussianity
	EG = zeros(ns, maxiter)

	# store error in the unmixing matrix -- only for testing (commented)
	# SSE = zeros(maxiter)


	# main loop
	t=0
	converged=false
	while !converged && t < maxiter
		t += 1
		copy!(Wp, W)
		deblend!(ica)

		if(ica.attrib == :complex)
			@simd for it=1:nt
				for is=1:ns
					ss=eps+abs(s[is, it])^2
					GWx[is, it] = log(ss)
					gWx[is, it] = inv(ss)
					dgWx[is, it] = -1. * inv(ss*ss)
				end
			end

			for ic = 1:ns
				b=zero(eltype(x))
				@simd for it=1:nt
					b+=gWx[ic,it]+(abs.(s[ic,it]).^2 * dgWx[ic,it])
				end
				b *= ntI
				for ic1=1:ns
					a=zero(eltype(x))
					@simd for it=1:nt
						a+=x[ic1,it]*conj(s[ic,it])*gWx[ic,it]
					end
					a *= ntI
					W[ic1, ic]=a-b*W[ic1,ic]
				end
			end
		elseif(ica.attrib == :real)
			@simd for it=1:nt
				for is=1:ns
					GWx[is,it]=exp(-0.5*s[is,it]*s[is,it])
					gWx[is,it]=s[is,it]*GWx[is,it]
					dgWx[is,it]=GWx[is,it]*(1.-(s[is,it]*s[is,it]))
				end
			end
			for ic = 1:ns
				b=zero(eltype(x))
				@simd for it=1:nt
					b+=dgWx[ic,it]
				end
				b *= ntI
				for ic1=1:ns
					a=zero(eltype(x))
					@simd for it=1:nt
						a+=x[ic1,it]*gWx[ic,it]
					end
					a *= ntI
					W[ic1, ic]=a-b*W[ic1,ic]
				end
			end
		else
			error("invalid ica.attrib")
		end


		# Symmetric decorrelation:
		copy!(W, W * sqrtm(inv(W'*W)))

		# compare W with Wp as a convergence criteria
		chgall=0.0
		for ic in 1:ns
			chg=0.0
			for ic1 in 1:ns
				chg += abs(W[ic1, ic]-Wp[ic1, ic]) 
			end
			chgall = max(chg, chgall) 
		end
		converged = (chgall < ica.tol)

		# storing the expectation of G(Wx) --- measure of Gaussianity
		for ic in 1:ns, it in 1:nt
			EG[ic,t] += (GWx[ic, it] * ntI)
		end

		# compute the distance between A and W'Q
		# this gives the measure of the success of the algorithm
		# note that W'Q and A can still be different by a 
		# permutation and scaling matrix, so...
		# if(!(A === nothing))
		#	SSE[t] = error_unmixing(W'*Q', A)
		# end
	end

	if(!(A === nothing))
	      SSE = error_unmixing(W'*Q', A)
	end
	
	# remove zeros
	EG = EG[:,1:t]
			
	# test the performance of the algorithm
	return ica, SSE, EG
end


"""
Work on it later, how to restart ICA if the starting model is bad.
"""
function fastica()
	tall=0
	converged_all=false

	while !converged_all && tall < 10
		tall += 1
		unmixing_initialize!(ica)
		converged_all = sum(abs.(EG[:,1]))/sum(abs.(EG[:,end]))
		println(converged_all)
	end
end




"""
Compute the error in the unmixing matrix.

* x is the unmixing matrix
* y is the original mixing matrix 

Ideally, x should be inverse of y upto scaling and permutation.
}	
"""
function error_unmixing{T}(x::Matrix{T}, y::Matrix{T}) 
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
Performs ICA for convolutive mixtures.
# Arguments
* `magic_recv`: a receiver index, where deblending is performed
* `recv_n`: total number of receivers
* `src_n`: total number of sources
* `grid`: `M1D` grid
* `nband`: number of frequency bins, where ICA is performed
* `nsubfac`: overlap factor (testing)
* `X`: blended data

# Output
* `Y`: deblended data at `magic_recv`
"""
function bica(;
	magic_recv::Int64=1,
	recv_n::Int64=1,
	src_n::Int64=1,
	grid::M1D=M1D(:samp1),
	nband::Int64=1, # number of bands 
	nsubfac::Int64=1, #
	X::Array{Float64}=zeros(grid.nx,recv_n) # data 
            	)

nband >= div(grid.nx,2) ? error("too many bands") :
nb = (mod(grid.nx, nband) == 0) ? div(grid.nx, nband) : error("change nsamples");

nbandtot = (nsubfac-1) * (nband-1) + nband
# initialize
Yband = zeros(src_n,nb)
Xband = zeros(recv_n,nb,nbandtot)
Mband = Array(ICA, nbandtot)

# divide into bands and do ICA in each band
for iband = 1: nbandtot
	# indices for bands
	ixmin = (iband-1) * div(nb,nsubfac) + 1
	ixmax = nb + (iband-1) * div(nb,nsubfac)

#	println("fffffff\t", ixmin, "\t", ixmax)

	# window out the data in the band
	Xband[:,:,iband] = transpose(X[ixmin:ixmax,:]);

	# do ICA
	Mband[iband] = ica(Xband[:,:,iband], src_n; 
	     do_whiten=true, fun=icagfun(:gaus), verbose=false, maxiter=200,
	     tol=1e-20, winit=ones(src_n,src_n))

	# calculate source signals temporarily
	Yband = transform(Mband[iband],Xband[:,:,iband]);

	# fix permutation problem
	Mband[iband].W = kurtosis(Yband[1,:]) < kurtosis(Yband[2,:]) ? 
			Mband[iband].W : Mband[iband].W[:,2:-1:1];

			#println("kurtosis 1\t", kurtosis(Yband[1,:]))
			#println("kurtosis 2\t", kurtosis(Yband[2,:]))
end

Yall = zeros(src_n,grid.nx);
Ybandprev = zeros(src_n,nb);
Yband = zeros(src_n,nb)
ixcmin  = 1
ixcmax  = nb - div(nb,nsubfac)
ixcumin = nb - div(nb,nsubfac) + 1
ixcumax = nb
ixcpmin = div(nb,nsubfac) + 1
ixcpmax = nb
for iband = 1: nbandtot 
	# indices of full band
	ixmin = (iband-1) * div(nb,nsubfac) + 1;
	ixmax = nb + (iband-1) * div(nb,nsubfac);

	# indices of the update zone that has no overlap
	ixminupdate = ixmax-div(nb,nsubfac)+1;
	ixmaxupdate = ixmax;


	Mbandtemp = Mband[iband];

	# minimum distortion principle -- implementation 1
#	WI = inv(Mbandtemp.W);
#	Mbandtemp.W *= diagm([WI[1,magic_recv],WI[2,magic_recv]]);


	# minimum distortion principle -- implementation 2
	W = transpose(Mbandtemp.W);
	A = inv(W);
	A[:,1] /= A[magic_recv,1];
	A[:,2] /= A[magic_recv,2];
	Mbandtemp.W = transpose(inv(A));



	#Mbandtemp.W[:,1] /= norm(Mbandtemp.W[:,1])
	#Mbandtemp.W[:,2] /= norm(Mbandtemp.W[:,2])

	# temporarily compute Yband 
	Yband = transform(Mbandtemp,Xband[:,:,iband]); 

	# normalize
#	Yband[1,:] /= norm(Yband[1,:])
#	Yband[2,:] /= norm(Yband[2,:])

	# update Yband 
	if((isequal(iband,1)) | (isequal(nsubfac,1)))
		Yall[1,ixmin:ixmax] = Yband[1,:];
		Yall[2,ixmin:ixmax] = Yband[2,:];

		# normalize
	#	Yall[1,ixmin:ixmax] /= norm(Yall[1,ixmin:ixmax])
#		Yall[2,ixmin:ixmax] /= norm(Yall[2,ixmin:ixmax])

	else
		# estimate scalars

		err1, sc1 = Misfits.error_after_scaling(Yband[1,ixcmin:ixcmax], Ybandprev[1,ixcpmin:ixcpmax])
		err2, sc2 = Misfits.error_after_scaling(Yband[2,ixcmin:ixcmax], Ybandprev[2,ixcpmin:ixcpmax])

		println("FFFF", err1, "\t",err2,"\t",sc1, "\t", sc2)
		# actual update
		Yall[1,ixminupdate:ixmaxupdate] = Yband[1,ixcumin:ixcumax] * sc1 ;
		Yall[2,ixminupdate:ixmaxupdate] = Yband[2,ixcumin:ixcumax] * sc2 ;

		# normalize after update
#		Yall[1,ixminupdate:ixmaxupdate] /= norm(Yall[1,ixminupdate:ixmaxupdate])
#		Yall[2,ixminupdate:ixmaxupdate] /= norm(Yall[2,ixminupdate:ixmaxupdate])

	end
	
	# save Yprev
	Ybandprev[1,:] = Yall[1,ixmin:ixmax];
	Ybandprev[2,:] = Yall[2,ixmin:ixmax];

end



#	# fix scaling problem using previous band
#	if((isequal(iband,1)) | (isequal(nsubfac,1)))
#		Mbandtemp = Mband[iband];
#
#		# minimum distortion principle 
#		WI = inv(Mbandtemp.W);
##		Mbandtemp.W *= diagm([WI[1,magic_recv],WI[2,magic_recv]]);
#if(!isequal(nsubfac,1))
#		Mbandtemp.W[:,1] /= norm(Mbandtemp.W[:,1])
#		Mbandtemp.W[:,2] /= norm(Mbandtemp.W[:,2])
#
#end
#
#		# final λbandtemp after fixing scaling and permutation
#		Yband = (transform(Mbandtemp,Xband[:,:,iband]));
#
#if(!isequal(nsubfac,1))
#	Yband[1,:] /= norm(Yband[1,:])
#	Yband[2,:] /= norm(Yband[2,:])
#end
#		Yall[:,ixmin:ixmax] = Yband[:,:];
#
#
#	else
#		λ = [dot(Ybandprev[1,ixcpmin:ixcpmax],Yband[1,ixcmin:ixcmax])/
#		     dot(Yband[1,ixcmin:ixcmax],    Yband[1,ixcmin:ixcmax]),
#		     dot(Ybandprev[2,ixcpmin:ixcpmax],Yband[2,ixcmin:ixcmax])/
#		     dot(Yband[2,ixcmin:ixcmax],    Yband[2,ixcmin:ixcmax])]
#
#		err1, sc1 =	Misfits.error_after_scaling(Yband[1,ixcmin:ixcmax], Ybandprev[1,ixcpmin:ixcpmax])
#		err2, sc2 =	Misfits.error_after_scaling( Yband[2,ixcmin:ixcmax], Ybandprev[2,ixcpmin:ixcpmax] )
#		println("GJJEHFFE", λ)
#		println(err1,"\t",sc1)
#		println(err2,"\t",sc2)
#
##		λ /= div(nb,nsubfac);
#		Mbandtemp = Mband[iband];
##		println(λ)
##		println("before corr", Mbandtemp.W, iband)
#		Mbandtemp.W = Mbandtemp.W * diagm(λ);
#		Mbandtemp.W[:,1] /= norm(Mbandtemp.W[:,1])
#		Mbandtemp.W[:,2] /= norm(Mbandtemp.W[:,2])
#		#Mbandtemp.W = Mbandtemp.W * diagm(λ./abs(λ));
#
#		#println("update in ",ixmax-div(nb,nsubfac)+1, "\t",ixmax)
#		#println("common ",ixcpmin,"\t", ixcpmax,"\t", ixcmin,"\t", ixcmax,"\t")
#
#		# final λbandtemp after fixing scaling and permutation
#		Yband = (transform(Mbandtemp,Xband[:,:,iband]));
#
#
#
#	Yall[1,ixmax-div(nb,nsubfac)+1: ixmax] = Yband[1,nb - div(nb,nsubfac) + 1:nb] * λ[1] ;
#	Yall[2,ixmax-div(nb,nsubfac)+1: ixmax] = Yband[2,nb - div(nb,nsubfac) + 1:nb] * λ[2] ;
#
#if(!isequal(nsubfac,1))
#	Yband[1,:] /= norm(Yband[1,:])
#	Yband[2,:] /= norm(Yband[2,:])
#end
#
#	end
#	println("after correc",Mbandtemp.W, iband)

#	Ybandprev = Yband;
#end




return reshape(Yall, (src_n,grid.nx))
end



"""
Convolutive mixing in the frequency domain.

# Arguments
* `As`:
* `Ab`:
* `B`:
* `S`:
* `fgrid`:

# Outputs
* `D`
* `d`
* `Ds`
* `ds`
* `Db`
* `db`
"""
function exact_freq_mixing(;
			   As::Array{Complex{Float64}}=nothing,
			   Ab::Array{Complex{Float64}}=nothing,
			   B::Array{Complex{Float64}}=nothing,
			   S::Array{Complex{Float64}}=nothing, 
			   fgrid::M1D=nothing
			  )
	recv_n = size(Ab,2)
	D = fill(complex(0.0,0.0),fgrid.nx,recv_n); 
	d = fill(complex(0.0,0.0),fgrid.nx,recv_n);
	Ds = fill(complex(0.0,0.0),fgrid.nx,recv_n); 
	ds = fill(complex(0.0,0.0),fgrid.nx,recv_n);
	Db = fill(complex(0.0,0.0),fgrid.nx,recv_n); 
	db = fill(complex(0.0,0.0),fgrid.nx,recv_n);
	for ir = 1:recv_n
		for iff = 1:fgrid.nx
			D[iff, ir] = Ab[iff, ir] * B[iff] + As[iff, ir] * S[iff]
			Db[iff, ir] = Ab[iff, ir] * B[iff] 
			Ds[iff, ir] = As[iff, ir] * S[iff]
		end
		d[:, ir] = ifft(D[:, ir])
		db[:, ir] = ifft(Db[:, ir])
		ds[:, ir] = ifft(Ds[:, ir])
	end
	return D, d, Ds, ds, Db, db
end



function do_whiten!(ica::ICA)
    nr, nt = ica.nr, ica.nt
    ntI = 1/nt
    x=ica.x
    xhat=ica.xhat
    Q=ica.Q
    C=Matrix{eltype(x)}(nr, nr)

    # construct covariance matrix
    A_mul_Bc!(C, x, x)
    scale!(C, ntI)

    # eigenvalue decomposition of the covariance matrix
    Efac = eigfact(C)
    ord = sortperm(Efac.values; rev=true)
    v, P = Efac.values[ord], Efac.vectors[:, ord]
    # whitening matrix
    copy!(Q, scale!(P, sqrt.(inv.(v))))

    # for testing --  this should be a diagonal matrix
    #println(W0' * C * W0)
    Ac_mul_B!(xhat, Q, x)
    return Q
end



end # module
