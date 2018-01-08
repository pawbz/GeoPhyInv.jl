using Revise
using JuMIT
using Base.Test
using ForwardDiff

ntgf = 5
nr = 10
nt = 25
gfobs=randn(ntgf, nr)
wavobs=randn(nt);
gfprecon, gfweights, wavprecon=JuMIT.DeConv.create_weights(ntgf,
		nt, gfobs, αexp=100., cflag=true)
pa=JuMIT.DeConv.Param(ntgf, nt, nr, gfobs=gfobs,
	wavobs=wavobs, wavnorm_flag=false, verbose=false,
	gfprecon=gfprecon, gfweights=gfweights, wavprecon=wavprecon,
	)

@time for attrib in [:gf, :wav]
	pa.attrib_inv=attrib
	x=randn(JuMIT.DeConv.ninv(pa))
	xa=similar(x)
	JuMIT.DeConv.x_to_model!(x, pa)
	JuMIT.DeConv.model_to_x!(xa, pa)
	# note that data is not linear when norm flag
	!pa.wavnorm_flag && @test x ≈ xa

	# compute Fx
	last_x=similar(x)
	JuMIT.DeConv.F!(pa, x)
	dcal=pa.cal.d

	# compute F^a dcala
	dcala=randn(size(pa.cal.d))
	gx=similar(x)
	JuMIT.DeConv.Fadj!(pa, x, gx, dcala)

	# note that data is not linear when norm_flag
	!pa.wavnorm_flag && @test dot(x, gx) ≈ dot(dcal, dcala)

	# compute gradient
	g1=similar(x);
	JuMIT.DeConv.func_grad!(g1, x,pa)
	g2=similar(x);
	f(x)=JuMIT.DeConv.func_grad!(nothing, x, last_x, pa)
	@time JuMIT.Inversion.finite_difference!(x -> JuMIT.DeConv.func_grad!(nothing, x, pa), x, g2, :central)
	@test g1 ≈ g2
end


# a simple DeConv test
ntgf = 50
nr = 40
nt = 5000

gfobs=randn(ntgf, nr)
wavobs=randn(nt)

pa=JuMIT.DeConv.Param(ntgf, nt, nr, gfobs=gfobs, wavobs=wavobs, wavnorm_flag=false,verbose=false)
storagewav=randn(size(pa.xwav))
storagegf=randn(size(pa.xgf))
# memory tests
using BenchmarkTools
println("===========")
pa.attrib_inv=:wav
@btime (randn!(pa.last_xwav); JuMIT.DeConv.F!(pa, pa.xwav,  ));
@btime JuMIT.DeConv.Fadj!(pa, pa.xwav, storagewav, pa.ddcal)
@btime JuMIT.DeConv.x_to_model!(pa.xwav, pa);
@btime JuMIT.DeConv.model_to_x!(pa.xwav, pa);
@btime JuMIT.DeConv.update_wav!(pa, pa.xwav,  pa.dfwav)
@btime JuMIT.DeConv.func_grad!(storagewav, pa.xwav,pa)
pba.attrib_inv=:gf
@btime JuMIT.DeConv.F!(pa, pa.xgf,);
@btime JuMIT.DeConv.Fadj!(pa, pa.xgf,storagegf, pa.ddcal)
@btime JuMIT.DeConv.x_to_model!(pa.xgf, pa);
@btime JuMIT.DeConv.model_to_x!(pa.xgf, pa);
@btime JuMIT.DeConv.update_gf!(pa, pa.xgf,  pa.dfgf)
@btime JuMIT.DeConv.func_grad!(storagegf, pa.xgf,  pa)


# check if there are any model discrepancies

begin
	gfobs=randn(ntgf, nr)
	wavobs=randn(nt)
	gfobs[1,:]=0.0; # mute
	gfprecon, gfweights, wavprecon=JuMIT.DeConv.create_weights(ntgf, nt, gfobs,
		αexp=10., cflag=true)
	paDeConv=JuMIT.DeConv.Param(ntgf, nt, nr, gfobs=gfobs, gfweights=gfweights,
	    fft_threads=false, wavnorm_flag=false, wavobs=wavobs, verbose=false, gfprecon=gfprecon,
	    wavprecon=wavprecon);
end

begin
	storagewav=randn(size(paDeConv.xwav))
	storagegf=randn(size(paDeConv.xgf))
	copy!(paDeConv.cal.wav, wavobs)
	copy!(paDeConv.cal.gf, gfobs)
end

begin
	paDeConv.attrib_inv=:wav
	JuMIT.DeConv.model_to_x!(paDeConv.xwav, paDeConv);
	@time JuMIT.DeConv.F!(paDeConv, paDeConv.xwav,  );
	@test paDeConv.cal.d ≈ paDeConv.obs.d
end

begin
	@time JuMIT.DeConv.func_grad!(storagewav, paDeConv.xwav, paDeConv)
	@test (storagewav ≈ zeros(storagewav))
end



begin
	paDeConv.attrib_inv=:gf
	JuMIT.DeConv.model_to_x!(paDeConv.xgf, paDeConv);
	@time JuMIT.DeConv.F!(paDeConv, paDeConv.xgf,  );
	@test paDeConv.cal.d ≈ paDeConv.obs.d
end
begin
	@time JuMIT.DeConv.func_grad!(storagegf, paDeConv.xgf, paDeConv)
	@test (storagegf ≈ zeros(storagegf))
end


# final test
@time JuMIT.DeConv.update_all!(pa)
f, α = JuMIT.Misfits.error_after_scaling(pa.cal.wav, wavobs)
@test (f<1e-3)
f, α = JuMIT.Misfits.error_after_scaling(pa.cal.gf, gfobs)
@test (f<1e-3)
