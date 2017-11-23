using Revise
using JuMIT
using Base.Test
using ForwardDiff

ntgf = 5
nr = 10
nt = 25
gfobs=randn(ntgf, nr)
gfobs[1,:]=0.0
wavobs=randn(nt);
gfprecon, gfweights, wavprecon=JuMIT.Decon.create_weights(ntgf,
		nt, gfobs, αexp=100., cflag=true)
pa=JuMIT.Decon.Param(ntgf, nt, nr, gfobs=gfobs,
	wavobs=wavobs, wavnorm_flag=true, verbose=false,
	gfprecon=gfprecon, gfweights=gfweights, wavprecon=wavprecon,
	)

@time for attrib in [:gf, :wav]
	pa.attrib_inv=attrib
	x=randn(JuMIT.Decon.ninv(pa))
	xa=similar(x)
	JuMIT.Decon.x_to_model!(x, pa)
	JuMIT.Decon.model_to_x!(xa, pa)
	# note that data is not linear when norm flag
	!pa.wavnorm_flag && @test x ≈ xa

	# compute Fx
	last_x=similar(x)
	JuMIT.Decon.F!(pa, x)
	dcal=pa.dcal

	# compute F^a dcala
	dcala=randn(size(pa.dcal))
	gx=similar(x)
	JuMIT.Decon.Fadj!(pa, x, gx, dcala)

	# note that data is not linear when norm_flag
	!pa.wavnorm_flag && @test dot(x, gx) ≈ dot(dcal, dcala)

	# compute gradient
	g1=similar(x);
	JuMIT.Decon.func_grad!(g1, x,pa)
	g2=similar(x);
	f(x)=JuMIT.Decon.func_grad!(nothing, x, last_x, pa)
	@time JuMIT.Inversion.finite_difference!(x -> JuMIT.Decon.func_grad!(nothing, x, pa), x, g2, :central)
	@test g1 ≈ g2
end


# a simple decon test
ntgf = 50
nr = 40
nt = 5000

gfobs=randn(ntgf, nr)
wavobs=randn(nt)

pa=JuMIT.Decon.Param(ntgf, nt, nr, gfobs=gfobs, wavobs=wavobs, wavnorm_flag=false,verbose=false)
storagewav=randn(size(pa.xwav))
storagegf=randn(size(pa.xgf))
# memory tests
for i in 1:4
	println("===========")
	pa.attrib_inv=:wav
	@time JuMIT.Decon.F!(pa, pa.xwav,  );
	@time JuMIT.Decon.Fadj!(pa, pa.xwav, storagewav, pa.ddcal)
	@time JuMIT.Decon.x_to_model!(pa.xwav, pa);
	@time JuMIT.Decon.model_to_x!(pa.xwav, pa);
	@time JuMIT.Decon.update_wav!(pa, pa.xwav, pa.last_xwav, pa.dfwav)
	@time JuMIT.Decon.func_grad!(storagewav, pa.xwav, pa.last_xwav, pa)
	pa.attrib_inv=:gf
	@time JuMIT.Decon.F!(pa, pa.xgf, pa.last_xgf, );
	@time JuMIT.Decon.Fadj!(pa, pa.xgf,storagegf, pa.ddcal)
	@time JuMIT.Decon.x_to_model!(pa.xgf, pa);
	@time JuMIT.Decon.model_to_x!(pa.xgf, pa);
	@time JuMIT.Decon.update_gf!(pa, pa.xgf, pa.last_xgf, pa.dfgf)
	@time JuMIT.Decon.func_grad!(storagegf, pa.xgf, pa.last_xgf, pa)
end


# check if there are any model discrepancies

begin
	gfobs=randn(ntgf, nr)
	wavobs=randn(nt)
	gfobs[1,:]=0.0; # mute
	gfprecon, gfweights, wavprecon=JuMIT.Decon.create_weights(ntgf, nt, gfobs, αexp=10., cflag=false)
	padecon=JuMIT.Decon.Param(ntgf, nt, nr, gfobs=gfobs, gfweights=gfweights,
	    fft_threads=false, wavnorm_flag=false, wavobs=wavobs, verbose=false, gfprecon=gfprecon,
	    wavprecon=wavprecon);
end

begin
	storagewav=randn(size(padecon.xwav))
	storagegf=randn(size(padecon.xgf))
	copy!(padecon.wav, wavobs)
	copy!(padecon.gf, gfobs)
end

begin
	padecon.attrib_inv=:wav
	JuMIT.Decon.model_to_x!(padecon.xwav, padecon);
	@time JuMIT.Decon.F!(padecon, padecon.xwav,  );
	@test padecon.dcal ≈ padecon.dobs
end

begin
	@time JuMIT.Decon.func_grad!(storagewav, padecon.xwav, padecon)
	@test (storagewav ≈ zeros(storagewav))
end




copy!(padecon.gf, gfobs)
padecon.attrib_inv=:gf
JuMIT.Decon.model_to_x!(padecon.xgf, padecon);
@time JuMIT.Decon.func_grad!(storagegf, padecon.xgf,  padecon)
@test (storagegf ≈ zeros(storagegf))


# final test
@time JuMIT.Decon.update_all!(pa)
f, α = JuMIT.Misfits.error_after_scaling(pa.wav, wavobs)
@test (f<1e-3)
f, α = JuMIT.Misfits.error_after_scaling(pa.gf, gfobs)
@test (f<1e-3)
