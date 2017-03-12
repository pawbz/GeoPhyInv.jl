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
nbh = Int64(floor(nb / 2));

# initialize
Yband = zeros(src_n,nb)
Xband = zeros(recv_n,nb,2*nband-1)
Mband = Array(ICA, 2*nband-1)

# divide into bands and do ICA in each band
for iband = 1: 2 * nband - 1
	# indices for bands
	ixmin = (div(iband,2))   * nb - (1-mod(iband,2)) * nbh + 1;
	ixmax = (div(iband,2)+1) * nb - (1-mod(iband,2)) * nbh;
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
	Mband[iband].W = entropy(Yband[1,:]) < entropy(Yband[2,:]) ? 
			Mband[iband].W : Mband[iband].W[:,2:-1:1];
end

Yall = zeros(src_n,grid.nx);
Ybandprev = zeros(src_n,nb);
Yband = zeros(src_n,nb)
for iband = 1: 2 * nband - 1 
	# indices
	ixmin = (div(iband,2))   * nb - (1-mod(iband,2)) * nbh + 1;
	ixmax = (div(iband,2)+1) * nb - (1-mod(iband,2)) * nbh;

	# temp 
	Yband = transform(Mband[iband],Xband[:,:,iband]); 

	# fix scaling problem using previous band
	if(iband > 1)
		λ = [sum(Ybandprev[1,nb-nbh:nb]./Yband[1,1:1+nbh])/nbh, 
			sum(Ybandprev[2,nb-nbh:nb]./Yband[2,1:1+nbh])/nbh]
		println(λ)
		Mbandtemp = Mband[iband];
		Mbandtemp.W *= (diagm(λ));
	#Mbandtemp.W *= diagm(diag((Mband[iband-1].W * inv(transpose(Mband[iband].W)))));
	else
		λ = [1., 1.]
		Mbandtemp = Mband[iband];

	end

	# final λbandtemp after fixing scaling and permutation
	Yband = (transform(Mbandtemp,Xband[:,:,iband]));
	Yall[:,ixmin:ixmax] = Yband;
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
