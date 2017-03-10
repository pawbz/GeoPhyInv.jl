module DSP

import SIT.Grid: M1D, M1D_npow2_tf
using Distributions

function get_random_tmax_signal(;
	fgrid::M1D=nothing, # frequency domain grid
	fmin::Float64=0.0, # minimum frequency
	fmax::Float64=nothing, # maximum frequency
	tmax::Float64=nothing, # frequency sampling is decided based on the length in time
	dist::Symbol=:gaussian # distribution type
	)

# initialize outputs
S = fill(complex(0.0,0.0),fgrid.nx);
s = fill(complex(0.0,0.0),fgrid.nx);

Δf = 1.0 / tmax
fvec = [f for f in fmin:Δf:fmax]
ifvec = fill(0, size(fvec))

for iff in eachindex(fvec)
	ifvec[iff] =  indmin((fgrid.x - fvec[iff]).^2.0)
end

if(dist == :gaussian)
	X = randn(size(fvec));
elseif(dist == :uniform)
	X = rand(Uniform(-1.0, 1.0), size(fvec))
else
	error("invalid dist")
end

sincvec = zeros(fgrid.x);
for iff in eachindex(fvec)
	sincvec = sinc(tmax.*(abs(fgrid.x - fgrid.x[ifvec[iff]])))
	S += complex(X[iff] .* sincvec, 0.0)
end

tgrid = M1D_npow2_tf(fgrid)

# remove mean (can as well put S(0) = 0)
s = ifft(S); s -=mean(s); S = fft(s);

return S, s




end



#function findfreq(; 
#		  x::Array{Float64},
#		  tgrid::Grid.Tim=Grid.Tim(),
#		  attrib::AbstractString="",
#		  threshold::Float64=1e-8
#		  )
#		  
#	dt, & ! time sampling interval
#	str, &  ! str  = PEAK returns peak frequency
#		! str  = MIN returns minumum frequency
#		! str  = MAX returns maximum frequency
#	threshold & ! threshold to select maximum and minumum frequency
#		    ! compared to the normalized power spectrum 
#)
#implicit none
#real, dimension(:), intent(in)          :: x
#real, intent(in)                        :: dt, threshold
#character(len=*), intent(in)		:: str
#optional				:: threshold
#
#real					:: thresholdin
#integer					:: npow2, nlag, ind,n
#real                                    :: findfreq
#real					:: df
#complex, dimension(:),allocatable       :: cx
#real, dimension(:), allocatable         :: ax, fvec
#logical, dimension(:), allocatable      :: mask
#
#if(present(threshold)) then; thresholdin = threshold; else; thresholdin = 1e-3; endif
#
#n=size(x); 
#call lenpow2(2*n,npow2)
#nlag = (npow2/2);
#df = 1.0 / npow2 / dt 
#allocate(cx(npow2), ax(npow2), fvec(npow2))
#cx = cmplx(rzero_de, rzero_de);
#ax = rzero_de;
#
#cx(1:n) = cmplx(x, rzero_de)
#
#call cfft(cx, npow2, 1)
#ax = real(cx * conjg(cx))
#fvec = npow2_vector(d=df, npow2=npow2)
#
#! normalizing
#if(maxval(ax) .eq. rzero_de) then
#        findfreq = 0.0
#        return
#else
#	ax = ax / maxval(ax)
#endif
#
#if(str .eq. 'PEAK') then
#	ind = maxloc(ax(1:nlag+1),dim=1)
#	findfreq = fvec(ind) 
#elseif(str .eq. 'MIN') then
#	allocate(mask(1:nlag+1))
#	mask = ax(1:nlag+1) .gt. thresholdin 
#	findfreq = minval(fvec(1:nlag+1), dim=1, mask=mask)
#	deallocate(mask)
#elseif(str .eq. 'MAX') then
#	allocate(mask(1:nlag+1))
#	mask = ax(1:nlag+1) .gt. thresholdin 
#	findfreq = maxval(fvec(1:nlag+1), dim=1, mask=mask)
#	deallocate(mask)
#else
#        findfreq = 0.0
#endif
#end function findfreq




end # module
