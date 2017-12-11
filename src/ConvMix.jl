__precompile__()


"""
Module to prepare convolutive mixtures.
And ICA is used for Blind Source Separation (BSS) of these mixtures.
"""
module ConvMix

import JuMIT.Conv
import JuMIT.CICA
import JuMIT.Grid
import JuMIT.Misfits
using DataFrames
using RecipesBase


type Param
	ns::Int # number of independent components
	obs::Vector{Conv.Param} # observed independent components
	nomix::Vector{Conv.Param} # unmixed data for simple xcorr
	mix::Vector{Conv.Param} # mixture, results of simple xcorr stored here
	"ica model"
	ica::CICA.ICA # results after performing ICA
	icaunmix::Vector{Conv.Param} # data after ica unmixing, this can be used for further processing
	"error in ge compared to gb"
	err::DataFrames.DataFrame
end


@userplot Plot


@recipe function f(p::Plot, rvec=nothing)
	pa=p.args[1]
end


"""
Greens functions are input to this constructor.
Source signal models are 
* `ntgf` : samples in Green function
* `nt` : samples in the long time series

* `ssparsepvec` : spartisty

"""
function Param(;ntgf=1, # instantaneos
	       nt=500,
	       nr=2,
	       ns=2,
	       gfobs=[randn(ntgf,nr) for is in 1:ns],
	       wavobs=[randn(nt) for is in 1:ns],
	       ica_func::Function=(x,ns)->CICA.ICA(x,ns),
	)

	obs=[Conv.Param(ntwav=nt, ntd=nt, ntgf=ntgf, dims=(nr,), wavlags=[nt-1, 0]) for is in 1:ns]
	nomix=deepcopy(obs); 
	mix=deepcopy(obs);
	icaunmix=deepcopy(obs);

	for is in 1:ns
		# copy observed Greens functions
		copy!(obs[is].gf,gfobs[is])
		# copy oroginal source wavelets
		for ir in 1:nr
			for i in 1:nt
				obs[is].wav[i,ir]=wavobs[is][i]
			end
		end
	end

	# model observed data for each source
	for is in 1:ns
		Conv.mod!(obs[is], :d)
	end


	# copy d and wav to nomix
	for is in 1:ns
		copy!(nomix[is].wav, obs[is].wav)
		copy!(nomix[is].d, obs[is].d)
	end

	# estimated gf without mixing and xcorr
	for is in 1:ns
		Conv.mod!(nomix[is], :gf)
	end

	#
	# copy sum of d and wav to mix
	for is in 1:ns
		for iss in 1:ns
			for i in eachindex(mix[is].d)
				mix[is].d[i] += obs[iss].d[i] 
			end
		end
		copy!(mix[is].wav, obs[is].wav)
	end

	# copy wavelets to icaunmix, replace this step with blind deconvolution?
	for is in 1:ns
		copy!(icaunmix[is].wav, obs[is].wav)
	end

	# estimated gf after mixing and xcorr
	for is in 1:ns
		Conv.mod!(mix[is], :gf)
	end


	x=transpose(rfft(mix[1].d, [1]))
	ica=ica_func(x, ns)

	err=DataFrame(rec_src=vec([(ir,is) for ir in 1:nr, is in 1:ns]), gf_no_mix_xcorr=zeros(nr*ns), gf_mix_xcorr=zeros(nr*ns), dat_mix_ica=zeros(nr*ns), 
	       gf_mix_ica=zeros(nr*ns))

	pa=Param(ns, obs, nomix,  mix, ica, icaunmix, err)

	update_err(pa)

	return pa

end

function update_err(pa)
	nr=size(pa.obs[1].gf,2)
	ns=pa.ns
	for is in 1:pa.ns
		for ir in 1:nr 
			d1=view(pa.nomix[is].gf, :, ir)
			d2=view(pa.mix[is].gf, :, ir)
			d0=view(pa.obs[is].gf, :, ir)
			pa.err[:no_mix_xcorr][ir+(is-1)*nr]=Misfits.error_after_scaling(d1, d0)[1]
			pa.err[:mix_xcorr][ir+(is-1)*nr]=Misfits.error_after_scaling(d2, d0)[1]
		end
	end
end
function Base.print(pa::Param)
	println("Convolutive Mixture\t")
	println("===================")
	println(pa.err)
end

function update_ica(cmix)
	x=rfft(cmix.d, [2])
	cmix.ica=ica_func(x, 2)
end


"""
Do deblending for every receiver individually.
And then update `icaunmix.d`
Then, `err_ica` is also updated with scaled least-squares error.
"""
function fastica!(pa::Param)
	nr=size(pa.obs[1].gf,2)
	ns=pa.ns

	for ir in 1:nr
		pa.ica.magic_recv=ir
		CICA.fastica!(pa.ica)

		s=transpose((pa.ica.s))
		s=irfft(s,pa.obs[1].ntd,[1])
		 
	end
	for is in 1:ns
		Conv.mod!(icaunmix[is], :gf)
	end
	for ir in 1:nr
		for is in 1:ns
			d=view(pa.icaunmix[is].d,:,ir)
			ss=view(s,:,is)
			copy!(d,ss)

			err=zeros(ns)
			for iss=1:ns
				d0=view(pa.obs[iss].d,:,ir)
				err[iss]=Misfits.error_after_scaling(d, d0)[1]
			end
			pa.err[:mix_ica][ir+(is-1)*nr]=minimum(err)
		end
	end
		

end





end
