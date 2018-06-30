module Smooth

using Conv
using Signals

#=
"""
! this subroutine returns a zero phase gaussian wi
ndow in time domain
n,& ! n must be either odd or a power of 2
! if n is odd, then the first nlag samples are negative lags
! followed by zero and positive lags
! if n is a power of 2, then lags are in wrap around order
! in this case one less negative lag is present, with gwin(1) as
! zero lag, followed by positive lags and negative lags are towards
! the end of the signal 
nwin, & ! this parameter controls the width of the gaussian window in
! signal, the width of an impulse after smoothin will be 2*nwin
gwin, & ! output zerophase gaussain window in time domain
gwin_freq& ! optional output in frequency domain [only when n = npow2]
  ! only real part is returned after cfft
)
"""
function gaussian_filter(nwin::Int64, nlag::Int64=2*nwin) 

	window = zeros(2*nlag+1)
	# time domain window
	if(nwin == 0)
		# if nwin = 0, then delta spike in time domain
		window[nlag + 1] =  1.0
	else
		# guassian window
		shrpcf = real(nwin)^(-2.0)
		for iw = 1: 2*nlag+1
			jsq          = ( nlag + 1 - iw )^2
			window[iw] = exp( -1.0 * shrpcf * real(jsq) )
		end
	end

	"""
	normalization of the filter such that the integral equals one
	necessary to preserve amplitudes after smoothing or whatever
	"""
	window ./= (sum(window))

end 

"""
! A zerophase gaussian window of length 2*nwin is designed
! and used to smooth datain
! width of an impulsive arrival in datout after smoothin will be 2*nwin
! The window itself can be seen as an the impluse response of the smoothing filter
! higher values of nwin means more smoothing
datin, & ! [n1*n2] input data matrix
datout, &! [n1*n2] output data after smoothing 
n1,&     ! first dimension of dat 
n2,&     ! second dimension of dat
nwin&    ! half of the width of the gaussian window used for smoothing
"""
function gaussian{N}(x::AbstractArray{Float64,N}, nwin::Vector{Int64})
	xsize=[size(x)...]

	pa=P_conv(dsize=xsize, ssize=ssize, gsize=xsize,)

	# convolution with a Gaussian window
	if(N==1)
		gwin=gaussian_filter(nwin[1])
		xin=deepcopy(x)
		Conv.conv!(x, xin, gwin, :d)
	elseif(N==2)
		gwin1=gaussian_filter(nwin[1])
		gwin2=gaussian_filter(nwin[2])
		xin=deepcopy(x)
		for i in 1:size(x,2)
			xx=view(x,:,i)
			xxin=view(xin,:,i)
			Conv.conv!(xx, xxin, gwin1, :d)
		end
		xin=deepcopy(x)
		for i in 1:size(x,1)
			xx=view(x,i,:)
			xxin=view(xin,i,:)
			Conv.conv!(xx, xxin, gwin2, :d)
		end
	else
		error("only implemented upto 2 dimensions")
	end
	return x
end


=#
function gaussian!{N}(xg::Array{Float64,N}, 
		      x::Array{Float64,N}, pa)
	Conv.mod!(pa, :d, d=xg, g=x)
end
end
