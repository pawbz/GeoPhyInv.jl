using JuMIT
using Base.Test

ntgf = 5
nr = 10
nt = 25
gfobs=randn(ntgf, nr)
wavobs=randn(nt);
pa=JuMIT.Decon.Param(ntgf, nt, nr, gfobs=gfobs, wavobs=wavobs, verbose=false)

@time for attrib in [:gf, :wav]
	pa.attrib_inv=attrib
	x=randn(JuMIT.Decon.ninv(pa))
	xa=similar(x)
	JuMIT.Decon.x_to_model!(x, pa)
	JuMIT.Decon.model_to_x!(xa, pa)
	@test x ≈ xa

	# compute Fx
	last_x=similar(x)
	JuMIT.Decon.F!(pa, x, last_x)
	dcal=pa.dcal

	# compute F^a dcala
	dcala=randn(size(pa.dcal))
	gx=similar(x)
	JuMIT.Decon.Fadj!(pa, gx, dcala)

	@test dot(x, gx) ≈ dot(dcal, dcala)
end


# a simple decon test
ntgf = 5
nr = 10
nt = 15

gfobs=randn(ntgf, nr)
wavobs=randn(nt)

pa=JuMIT.Decon.Param(ntgf, nt, nr, gfobs=gfobs, wavobs=wavobs, verbose=false)
storagewav=randn(size(pa.xwav))
storagegf=randn(size(pa.xgf))
# memory tests
for i in 1:4
	println("===========")
	pa.attrib_inv=:wav
	@time JuMIT.Decon.F!(pa, pa.xwav, pa.last_xwav, );
	@time JuMIT.Decon.Fadj!(pa, storagewav, pa.ddcal)
	@time JuMIT.Decon.x_to_model!(pa.xwav, pa);
	@time JuMIT.Decon.model_to_x!(pa.xwav, pa);
	@time JuMIT.Decon.update_wav!(pa, pa.xwav, pa.last_xwav, pa.dfwav)
	@time JuMIT.Decon.func_grad!(storagewav, pa.xwav, pa.last_xwav, pa)
	pa.attrib_inv=:gf
	@time JuMIT.Decon.F!(pa, pa.xgf, pa.last_xgf, );
	@time JuMIT.Decon.Fadj!(pa, storagegf, pa.ddcal)
	@time JuMIT.Decon.x_to_model!(pa.xgf, pa);
	@time JuMIT.Decon.model_to_x!(pa.xgf, pa);
	@time JuMIT.Decon.update_gf!(pa, pa.xgf, pa.last_xgf, pa.dfgf)
	@time JuMIT.Decon.update_wav!(pa, pa.xwav, pa.last_xwav, pa.dfwav)
	@time JuMIT.Decon.func_grad!(storagegf, pa.xgf, pa.last_xgf, pa)
end


# final test
@time JuMIT.Decon.update_all!(pa)
f, α = JuMIT.Misfits.error_after_scaling(pa.wav, wavobs)
@test (f<1e-3)
f, α = JuMIT.Misfits.error_after_scaling(pa.gf, gfobs)
@test (f<1e-3)
