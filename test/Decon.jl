using Revise
using JuMIT
using Base.Test
using ForwardDiff

ntgf = 5
nr = 10
nt = 25
gfobs=randn(ntgf, nr)
wavobs=randn(nt);
pa=JuMIT.Decon.Param(ntgf, nt, nr, gfobs=gfobs, wavobs=wavobs, wavnorm_flag=true, verbose=false)

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
	JuMIT.Decon.F!(pa, x, last_x)
	dcal=pa.dcal

	# compute F^a dcala
	dcala=randn(size(pa.dcal))
	gx=similar(x)
	JuMIT.Decon.Fadj!(pa, x, gx, dcala)

	# note that data is not linear when norm_flag
	!pa.wavnorm_flag && @test dot(x, gx) ≈ dot(dcal, dcala)

	# compute gradient
	g1=similar(x);
	JuMIT.Decon.func_grad!(g1, x, last_x,pa)
	g2=similar(x);
	f(x)=JuMIT.Decon.func_grad!(nothing, x, last_x, pa)
	@time JuMIT.FWI.finite_difference!(x -> JuMIT.Decon.func_grad!(nothing, x, last_x, pa), x, g2, :central)
	@test g1 ≈ g2
end


# a simple decon test
ntgf = 150
nr = 40
nt = 15000

gfobs=randn(ntgf, nr)
wavobs=randn(nt)

pa=JuMIT.Decon.Param(ntgf, nt, nr, gfobs=gfobs, wavobs=wavobs, wavnorm_flag=false,verbose=false)
storagewav=randn(size(pa.xwav))
storagegf=randn(size(pa.xgf))
# memory tests
for i in 1:4
	println("===========")
	pa.attrib_inv=:wav
	@time JuMIT.Decon.F!(pa, pa.xwav, pa.last_xwav, );
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


# final test
@time JuMIT.Decon.update_all!(pa)
f, α = JuMIT.Misfits.error_after_scaling(pa.wav, wavobs)
@test (f<1e-3)
f, α = JuMIT.Misfits.error_after_scaling(pa.gf, gfobs)
@test (f<1e-3)
