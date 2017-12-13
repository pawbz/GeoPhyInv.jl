# blind deconvolution
__precompile__()

module DeConvOP
import JuMIT.DSP
import JuMIT.Operators
import JuMIT.Inversion
import JuMIT.Misfits
import JuMIT.Inversion
using Optim, LineSearches
using RecipesBase
using StatsBase

#=
	ntgf	length of greens
	nt	length of wav and data
	np2	zero padded length
	g	[ntgf*nr] greens functions
	Pg	[np2*nr, ntgf*nr] padding matrix greens
	Pw	[np2*nr, nt] padding matrix wav

	=#

type Param
	nt::Int64
	nr::Int64
	ntgf::Int64
	gfobs::Array{Float64,1}
	gf::Array{Float64,1}
	gfprecon::Array{Float64,1}
	gfweights::Array{Float64,1}
	dobs::Array{Complex{Float64},1}
	dcal::Array{Complex{Float64},1}
	wavobs::Vector{Float64}
	wav::Vector{Float64}
	wavprecon::Vector{Float64}
	G::Array{Complex{Float64},1}
	W::Array{Complex{Float64},1}
	Dw::Array{Complex{Float64},2}
	Dg::Array{Complex{Float64},2}
	DPg::Matrix{Complex{Float64}}
	DPw::Matrix{Complex{Float64}}
	iDT::Matrix{Complex{Float64}}
	SW::SparseMatrixCSC{Complex{Float64},Int64}
	SG::SparseMatrixCSC{Complex{Float64},Int64}
	Og::Array{Complex{Float64},2}
	Ow::Array{Complex{Float64},2}
	attrib_inv::Symbol
	verbose::Bool
end

@userplot Plot


@recipe function f(p::Plot, rvec=nothing)
	pa=p.args[1]
	(rvec===nothing) && (rvec=1:pa.nr)

	# time vectors
	# autocorr wav
	awavobs=autocor(pa.wavobs, 1:pa.nt-1, demean=true)
	awav=autocor(pa.wav, 1:pa.nt-1, demean=true)
	wavli=max(maximum(abs,awavobs), maximum(abs,awav))
	# autocorr gf 
	agfobs=autocor(pa.gfobs,1:pa.ntgf-1, demean=true)
	agf=autocor(pa.gf,1:pa.ntgf-1, demean=true)
	gfli=max(maximum(abs,agfobs), maximum(abs,agf))

	# cut receivers
	dcal=pa.dcal[:,rvec]
	dobs=pa.dobs[:,rvec]

	layout := (5,3)

	@series begin        
		subplot := 1
#		aspect_ratio := :auto
		legend := false
		pa.wavobs
	end
	@series begin        
		subplot := 2
		legend := false
		pa.gfobs
	end
#
#	@series begin        
#		subplot := 3
#		legend := false
#		pa.dobs
#	end
#
#	@series begin        
#		subplot := 4
#		legend := false
#		pa.wav
#	end
#	@series begin        
#		subplot := 5
#		legend := false
#		pa.gf[:,rvec]
#	end
#
#	@series begin        
#		subplot := 6
#		legend := false
#		pa.dcal[:,rvec]
#	end
#
#	@series begin        
#		subplot := 7
#		legend := false
#		awavobs
#
#		
#	end
#	@series begin        
#		subplot := 8
#		legend := false
#		agfobs[:, rvec]
#	end
#
#	@series begin        
#		subplot := 9
#		legend := false
#		
#	end
#
#	@series begin        
#		subplot := 10
#		legend := false
#		awav
#
#		
#	end
#	@series begin        
#		subplot := 11
#		legend := false
#		agf[:,rvec]
#	end
#
#	@series begin        
#		subplot := 12
#		legend := false
#		pa.dcal[:,rvec]-pa.dobs[:,rvec]
#	end
#
#
#
#
#	@series begin        
#		subplot := 13
#		aspect_ratio := :equal
#		seriestype := :histogram2d
#		title := "Scatter Wav"
#		legend := false
#		awavobs, awav
#	end
#
#	@series begin        
#		subplot := 14
#		aspect_ratio := :equal
#		seriestype := :histogram2d
#		title := "Scatter Wav"
#		legend := false
#		agfobs, agf
#	end
#
#	@series begin        
#		subplot := 15
#		aspect_ratio := 1
#		seriestype := :histogram2d
#		title := "Scatter Wav"
#		legend := false
#		pa.dobs, pa.dcal
#	end
#


end

"""
`gfprecon` : a preconditioner applied to each Greens functions [ntgf]
"""
function Param(ntgf, nt, nr; 
	       gfprecon=nothing,
	       gfweights=nothing,
	       gfoptim=nothing,
	       gfαvec=nothing,
	       wavoptim=nothing,
	       wavαvec=nothing,
	       wavprecon=nothing,
	       wavnorm_flag=false,
	       fft_threads=false,
	       dobs=nothing, gfobs=nothing, wavobs=nothing, verbose=false, attrib_inv=:gf,) 
	(dobs===nothing) && (dobs=complex.(zeros(nt*nr)))
	dcal=complex.(zeros(nt*nr))
	

		# initial values are random
	wav=randn(nt)
	gf=randn(ntgf*nr)

	# fft dimension for plan
	np2=maximum([2*nt, 2*ntgf])

	# padding array for 
	Pg=Operators.splitpad(Complex128, nr*ntgf, np=np2-ntgf, nfrac=nr)
	Pd=Operators.splittrun(Complex128, nr*nt, np=np2-nt, nfrac=nr)
	Pw=Operators.pad(Complex128, nt, np=np2-nt, ntimes=nr)

	DFT=fft(eye(np2), 1)
	DPg=Operators.blkdiag(Complex128,DFT,nr)*Pg
	DPw=Operators.blkdiag(Complex128,DFT,nr)*Pw

	
	#DPw=kron(speye(nr),fft(speye(np2), 1))*Pw # perfrom fft and zero padding
	iDFT=ifft(eye(np2),1)
	iDT=Pd*Operators.blkdiag(Complex128,iDFT, nr)

	#iDT=Pd*kron(eye(nr),ifft(eye(np2), 1)) # perform ifft and truncation
	
#	Q=spdiagm(complex.(zeros(np2))) # allocate diagonal operator

	G=complex.(zeros(np2*nr))
	Dw=complex.(zeros(np2*nr, ntgf*nr))
	Dg=complex.(zeros(np2*nr, nt))
	W=complex.(zeros(np2*nr))

	# sparse diagonal matrix for wav in frequency domain
	SW=spzeros(Complex{Float64},size(W,1), size(W,1));
	SG=spzeros(Complex{Float64},size(G,1), size(G,1));

	Ow=complex.(zeros(nt*nr, nr*ntgf)) # real operator that has to be applied on gf to get d
	
        Og=complex.(zeros(nt*nr, nt)) # real operator that has to be applied on wav to get d

	# create dummy gfobs if necessary
	(gfobs===nothing) && (gfobs=zeros(gf))
	# create dummy wavobs if necessary
	(wavobs===nothing) && (wavobs=zeros(wav))

	# create gf precon
	(gfprecon===nothing) && (gfprecon=ones(gf))
	(gfweights===nothing) && (gfweights=ones(gf))

	# create gf precon
	(wavprecon===nothing) && (wavprecon=ones(wav))

	pa = Param(nt, nr, ntgf, gfobs, gf, gfprecon, gfweights, dobs, dcal, wavobs, wav, wavprecon, G, W, Dw, Dg,
	    DPg, DPw, iDT, SW, SG, Og, Ow,
	    attrib_inv, verbose)
 
	if iszero(pa.dobs) 
		((gfobs===nothing) | (wavobs===nothing)) && error("need gfobs and wavobs")
		F!(pa,gf=gfobs, wav=wavobs,d=pa.dobs) 
	end

	initialize!(pa)

	return pa
	
end

"""
Green functions after zero padding and FFT
"""
function update_G!(pa, gf)
	A_mul_B!(pa.G, pa.DPg, gf)
end

"""
Wavelet after zero padding and FFT
"""
function update_W!(pa, wav)
	A_mul_B!(pa.W, pa.DPw, wav)
end


"""
Prepare sparse matrix for Wavelet
"""
function update_SW!(pa, wav)
	update_W!(pa, wav)
	for i in eachindex(pa.W)
		pa.SW[i,i]=pa.W[i]
	end
end

"""
Prepare sparse matrix for Greens functions
"""
function update_SG!(pa, gf)
	update_G!(pa, gf)
	for i in eachindex(pa.G)
		pa.SG[i,i]=pa.G[i]
	end
end

"""
Operator to be applied on gf
"""
function update_Ow!(pa, wav)
	update_SW!(pa, wav)
	A_mul_B!(pa.Dw, pa.SW, pa.DPg)
	A_mul_B!(pa.Ow,pa.iDT,pa.Dw)
end


function update_Og!(pa, gf)
	update_SG!(pa, gf)
	A_mul_B!(pa.Dg, pa.SG, pa.DPw)
	A_mul_B!(pa.Og,pa.iDT,pa.Dg)
end

function F!(pa;gf=pa.gf,wav=pa.wav,d=pa.dcal)
	update_Og!(pa,gf)
	A_mul_B!(d, pa.Og, wav)
end

function ninv(pa)
	if(pa.attrib_inv == :wav)
		return pa.nt
	else(pa.attrib_inv == :gf)
		return pa.ntgf*pa.nr
	end
end


function error(pa) 
	fwav = Misfits.error_after_normalized_autocor(pa.wav, pa.wavobs)
	fgf = Misfits.error_after_normalized_autocor(pa.gf, pa.gfobs)
	f = Misfits.error_squared_euclidean!(nothing, real.(pa.dcal), real.(pa.dobs), nothing, norm_flag=true)

	println("Blind Decon\t")
	println("===========")
	println("error in estimated wavelet:\t", fwav)
	println("error after autocor in estimated Green Functions:\t", fgf)
	println("normalized error in the data:\t", f)

	return fwav, fgf, f
end 


function update_gf!(pa)
	update_Ow!(pa,pa.wav)
	gf=real.(pinv(pa.Ow)*pa.dobs)
	copy!(pa.gf,gf)
	F!(pa)

	return Misfits.error_squared_euclidean!(nothing, real.(pa.dcal), real.(pa.dobs), nothing, norm_flag=true)
end

function update_wav!(pa)
	update_Og!(pa,pa.gf)
	wav=real.(pinv(pa.Og)*pa.dobs)
	copy!(pa.wav, wav)
	F!(pa)
	return Misfits.error_squared_euclidean!(nothing, real.(pa.dcal), real.(pa.dobs), nothing, norm_flag=true)

end

function remove_gfprecon!(pa)
	for i in eachindex(pa.gfprecon)
		if(pa.gfprecon[i]≠0.0)
			pa.gfprecon[i]=1.0
		end
	end
end

"""
* re_init_flag :: re-initialize inversions with random input or not?
"""
function update_all!(pa; max_roundtrips=100, max_reroundtrips=10, ParamAM_func=nothing, roundtrip_tol=1e-3)

	if(ParamAM_func===nothing)
		ParamAM_func=x->Inversion.ParamAM(x, optim_tols=[1e-5, 1e-5],name="Blind Decon",
				    roundtrip_tol=roundtrip_tol, max_roundtrips=max_roundtrips,
				    max_reroundtrips=max_reroundtrips,
				    min_roundtrips=10,
				    reinit_func=x->initialize!(pa))
	end

	
	# create alternating minimization parameters
	f1=x->update_wav!(pa)
	f2=x->update_gf!(pa)
	paam=ParamAM_func([f1, f2])

	# do inversion
	Inversion.go(paam)

	# print errors
	DeconOP.error(pa)
end


function initialize!(pa)
	# starting random models
	randn!(pa.wav)
	randn!(pa.gf)
end





"""
Create preconditioners using the observed Green Functions.
* `cflag` : causaulity flag
"""
function create_weights(ntgf, nt, gfobs; αexp=0.0, cflag=true)

	nr=size(gfobs,2)
	wavprecon=ones(nt)
	gfprecon=ones(ntgf, nr); 
	gfweights=ones(ntgf, nr); 
	minindz=ntgf
	gfweights=ones(ntgf, nr)
	for ir in 1:nr
		gf=normalize(view(gfobs,:,ir))
		indz=findfirst(x->abs(x)>1e-6, gf)
		if(cflag)
			if(indz==0)
				gfprecon[:,ir]=0.0    
				gfweights[:,ir]=0.0
			else
				gfprecon[1:indz-1,ir]=0.0    
				gfweights[1:indz-1,ir]=0.0
#				wavprecon[end-indz+1:end]=0.0
			end
		end
		if(indz≠0)
			for i in indz+1:ntgf        
				gfweights[i,ir]=exp(αexp*(i-indz-1)/ntgf)  # exponential weights
				gfprecon[i,ir]=exp(αexp*(i-indz-1)/ntgf)  # exponential weights
			end
		end
	end

	return gfprecon, gfweights, wavprecon

end

end # module
