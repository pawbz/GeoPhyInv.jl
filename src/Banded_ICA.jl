module Banded_ICA

using MultivariateStats
using Distributions
import SIT.Grid:M1D

function bica(;
	recv_n::Int64=1,
	src_n::Int64=1,
	grid::M1D=M1D(:samp1),
	nband::Int64=1, # number of bands 
	X::Array{Float64}=zeros(grid.nx,recv_n) # data 
            	)

nb = (mod(grid.nx, nband) == 0) ? div(grid.nx, nband) : error("change nsamples");

nsubfac = 8;
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

	println("fffffff\t", ixmin, "\t", ixmax)

	# window out the data in the band
	Xband[:,:,iband] = transpose(X[ixmin:ixmax,:]);

	# do ICA
	Mband[iband] = fit(ICA, Xband[:,:,iband], src_n; 
	     do_whiten=true, fun=icagfun(:gaus), verbose=false, maxiter=50,
	     tol=1e-20, winit=ones(src_n,src_n))

	# calculate source signals temporarily
	Yband = transform(Mband[iband],Xband[:,:,iband]);

	# fix permutation problem
	Mband[iband].W = kurtosis(Yband[1,:]) < kurtosis(Yband[2,:]) ? 
			Mband[iband].W : Mband[iband].W[:,2:-1:1];

			println("kurtosis 1\t", kurtosis(Yband[1,:]))
			println("kurtosis 2\t", kurtosis(Yband[2,:]))
end

Yall = zeros(src_n,grid.nx);
Ybandprev = zeros(src_n,nb);
Yband = zeros(src_n,nb)
ixcmin  = 1
ixcmax  = nb - div(nb,nsubfac)
ixcpmin = div(nb,nsubfac) + 1
ixcpmax = nb
for iband = 1: nbandtot 
	# indices
	ixmin = (iband-1) * div(nb,nsubfac) + 1
	ixmax = nb + (iband-1) * div(nb,nsubfac)

	# temp 
	Yband = transform(Mband[iband],Xband[:,:,iband]); 

	# fix scaling problem using previous band
	if(isequal(iband,1))
		Mbandtemp = Mband[iband];
		# final λbandtemp after fixing scaling and permutation
		Yband = (transform(Mbandtemp,Xband[:,:,iband]));
		Yall[:,ixmin:ixmax] = Yband[:,:];
	else
		λ = [dot(Ybandprev[1,ixcpmin:ixcpmax],Yband[1,ixcmin:ixcmax])/
		     dot(Yband[1,ixcmin:ixcmax],    Yband[1,ixcmin:ixcmax]),
		     dot(Ybandprev[2,ixcpmin:ixcpmax],Yband[2,ixcmin:ixcmax])/
		     dot(Yband[2,ixcmin:ixcmax],    Yband[2,ixcmin:ixcmax])]
		Mbandtemp = Mband[iband];
		println(λ)
		println("before corr", Mbandtemp.W, iband)
		#Mbandtemp.W = Mbandtemp.W * diagm(λ);
		Mbandtemp.W = Mbandtemp.W * diagm(λ./abs(λ));
		#Mbandtemp.W *= diagm(diag((Mband[iband-1].W * inv(transpose(Mband[iband].W)))));

		# final λbandtemp after fixing scaling and permutation
		Yband = (transform(Mbandtemp,Xband[:,:,iband]));
		Yall[:,ixmax-div(nb,nsubfac)+1: ixmax] = Yband[:,nb - div(nb,nsubfac) + 1:nb];
	end
	println("after correc",Mbandtemp.W, iband)

	Ybandprev = Yband;
end




return reshape(Yall, (src_n,grid.nx))
end

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
	for ir = 1:recv_n
		for iff = 1:fgrid.nx
			D[iff, ir] = Ab[iff, ir] * B[iff] + As[iff, ir] * S[iff]
		end
		d[:, ir] = ifft(D[:, ir])
	end
	return D, d
end


end # module
