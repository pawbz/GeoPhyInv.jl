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
@time JuMIT.Decon.update_all!(pa)
f, α = JuMIT.Misfits.error_after_scaling(pa.wav, wavobs)
@test (f<1e-3)
f, α = JuMIT.Misfits.error_after_scaling(pa.gf, gfobs)
@test (f<1e-3)
