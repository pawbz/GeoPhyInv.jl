module Banded_ICA

using MultivariateStats
import SIT.Grid
using Distributions

function bica(;
	recv_n::Int64=1,
	src_n::Int64=1,
	grid::Grid.M1D=Grid.M1D(:samp1),
	nband::Int64=1, # number of bands 
	X::Array{Float64}=zeros(grid.nx,recv_n) # data 
            	)

nb = (mod(grid.nx, nband) == 0) ? div(grid.nx, nband) : error("change nsamples")

# initialize
Yall = zeros(src_n,nb,nband);
Yband = zeros(src_n,nb)
Qref = zeros(src_n, recv_n)
# divide into bands and do ICA in each band
for iband = 1:nband
	# window out the data in the band
	Xband = transpose(X[(iband-1)*nb+1:iband*nb,:]);
	Mband = fit(ICA, Xband, src_n; 
	     do_whiten=true, fun=icagfun(:gaus), verbose=false, maxiter=50,
	     tol=1e-20, winit=ones(src_n,src_n))

	# calculate source signals temporarily
	Yband = transform(Mband,Xband);

	# fix permutation problem
	Mband.W = kurtosis(Yband[1,:]) < kurtosis(Yband[2,:]) ? Mband.W : Mband.W[:,2:-1:1];
	
	# save W during the first band as reference Q
	if(isequal(iband,1)); Qref = transpose(Mband.W); end;

	# fix scaling problem using Qref
	Mband.W = Mband.W * diagm(diag((Qref * inv(transpose(Mband.W)))));

	# final λbandtemp after fixing scaling and permutation
	Yall[:,:,iband] = (transform(Mband,Xband));
end

return reshape(Yall, (src_n,grid.nx))
end

end # module
