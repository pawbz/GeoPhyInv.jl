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
	"blended data (nrec, nt)"
	x::Array{T,2}
	"number of sources or independent components"
	ns::Int64
	"store G(Wx)"
	GWx::Matrix{T}
	"store g(Wx)"
	gWx::Matrix{T}
	"store dg(Wx)"
	dgWx::Matrix{T}
	"unmixing matrix --- this is updated after ICA"
	W::Array{T}
	"previous unmixing matrix"
	Wp::Array{T}
	"store mean vector"
	mv::Vector{T}
	"real or complex ica"
	attrib::Symbol
	"maximum iterations"
	maxiter:Int64
	"tolerance"
	tol::Float64
end

function ICA(x,
	     maxiter::Int=100,
	     tol::Real=1.0e-10,                 # convergence tolerance
	     Win=nothing,  # init guess of W, size (m, k)
	     )

	T = eltype(x)

	# check input arguments
	nrec, nt = size(x)
	(nt < 1) && error("there must be more than one samples, i.e. nt > 1.")
	(nsrc > nrec) && error("nsrc must not exceed nrec")

	# allocating vectors, that store G, g, and dg
	GWx = Matrix{T}(nc, nt)
	gWx = Matrix{T}(nc, nt)
	dgWx = Matrix{T}(nc, nt)

	# random initialization of unmixing matrix
	if(T == Float64)
		attrib=:real
		W = (Win === nothing) ? randn(nc,nc) : Win
	elseif(T == Complex{Float64})
		attrib=:complex
		W = (Win === nothing) ? complex.(randn(nc,nc) + randn(nc,nc)) : Win
	else
	    	error("unknown type")
	end
	Wp = similar(W)


	return ICA(x, )

end

# Fast ICA
#
# Reference:
#
#   Aapo Hyvarinen and Erkki Oja
#   Independent Component Analysis: Algorithms and Applications.
#   Neural Network 13(4-5), 2000.
#
function fastica!{T}(x::Matrix{T}, 
		     nsrc::Int;
		     whiten_flag::Bool=true,             # whether to perform pre-whitening 
                     tol::Real=1.0e-10,                 # convergence tolerance
		     A=nothing, # mixing matrix used to estimate performance of algorithm (unknown)
                          verbose::Bool=false)              # whether to display iterations

	(maxiter < 1) && error("maxiter must be greater than 1.")
	(tol < 0.) && error("tol must be positive.")

	# preprocess data
	mv = Vector{T}(nrec)

	# calculate and store mean
	for irec=1:nrec
		mv[irec] =  mean(x[irec, :])
	end
	# remove mean from x; will be added back later
	for irec=1:nrec
		copy!(x[irec, :], (x[irec, :] - mv[irec]))
	end
	if(whiten_flag)
		Q = do_whiten!(x)
	end

	nc, nt = size(x)
	# normalize each column
	for ic = 1:nc
		w = view(W,:,ic)
		scale!(w, 1.0 / sqrt(sum(abs2, w)))
	end

	# for complex ICA
	eps = 0.1

	# preallocating GWx, gWx, dgWx
	GWx = Matrix{T}(nc, nt)
	gWx = Matrix{T}(nc, nt)
	dgWx = Matrix{T}(nc, nt)

	EG = zeros(nc, maxiter)
	SSE = zeros(maxiter)
	# main loop
	t = 0
	converged = false
	while !converged && t < maxiter
		t += 1
		Wp = copy(W)
	    	for ic = 1:nc
			if(attrib_ica == :complex)
				xx = W[:,ic]'*x
				y = abs.(xx).^2
				attrib_func = :log
				if(attrib_func == :log)
					GWx[ic,:] = log.(eps + y)
					gWx[ic,:] = 1./(eps + y)
					dgWx[ic,:] = -1./(eps + y).^2
				elseif(attrib_func == :gaus)
					GWx[ic,:] = -exp.(-0.5*xx.*xx)
					gWx[ic,:] = xx.*exp.(-0.5*xx.*xx)
					dgWx[ic,:] = exp.(-0.5*xx.*xx) - xx.*xx
				else
					error("invalid attrib_func")
				end

				a = mean(x .* (ones(nc,1)*conj(W[:,ic]'*x)) .* (ones(nc,1)*gWx[[ic],:]),2)[:,1]
				b = mean(gWx[[ic],:] + abs.(W[:,ic]'*x).^2 .* dgWx[[ic],:]) * W[:,ic]
				W[:,ic] = a - b

			elseif(attrib_ica == :real)
				xx = W[:,ic]'*x
				GWx[ic,:] = -exp.(-0.5*xx.*xx)
				gWx[ic,:] = xx.*exp.(-0.5*xx.*xx)
				dgWx[ic,:] = exp.(-0.5*xx.*xx) - xx.*xx
				a = mean(x .* (ones(nc,1)*gWx[[ic],:]),2)[:,1]
				b = mean(dgWx[[ic],:]) * W[:,ic]
				W[:,ic] = a - b
			else
				error("invalid attrib_ica")
			end

		end

		# Symmetric decorrelation:
		copy!(W, W * sqrtm(inv(W'*W)))

		# compare W with Wp as a convergence criteria
		converged = (maximum(abs.(W-Wp)) < tol)

		# storing the expectation of G(Wx) --- measure of Gaussianity
		EG[:,t] = mean(GWx,2);

		# compute the distance between A and W'Q
		# this gives the measure of the success of the algorithm
		# note that W'Q and A can still be different by a 
		# permutation and scaling matrix, so...
		if(!(A === nothing))
			SSE[t] = Misfits.error_after_scaling_permutation(W'*Q, A)
		end
	end

	# demixing 
	x = W' * x

	# construct model
	if(whiten_flag)
		W = W' * Q
	end

	# adding mean back to x; avoided making a copy of x
	for irec=1:nrec
		copy!(x[irec, :], (x[irec, :] + mv[irec]))
	end

	return x, W, EG[:,1:t], SSE[1:t]
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



# Independent Component Analysis

centralize(x::AbstractVector, m::AbstractVector) = (isempty(m) ? x : x - m)::typeof(x)
centralize(x::AbstractMatrix, m::AbstractVector) = (isempty(m) ? x : x .- m)::typeof(x)

decentralize(x::AbstractVector, m::AbstractVector) = (isempty(m) ? x : x + m)::typeof(x)
decentralize(x::AbstractMatrix, m::AbstractVector) = (isempty(m) ? x : x .+ m)::typeof(x)


preprocess_mean{T<:Number}(X::Matrix{T}, m) = (m == nothing ? vec(Base.mean(X, 2)) :
                                                      m == 0 ? T[] : 
                                                      m)::Vector{T}
function do_whiten!(x)
    nrec, nt = size(x)
    # construct covariance matrix
    C = (x * x') ./ nt
    # eigenvalue decomposition of the covariance matrix
    Efac = eigfact(C)
    v, P = Efac.values, Efac.vectors
    # whitening matrix
    W0 = scale!(P, sqrt.(inv.(v)))
    # for testing --  this should be a diagonal matrix
    #println(W0' * C * W0)
    copy!(x, W0' * x)
    return W0'
end



end # module
