using Statistics
using Test
using GeoPhyInv



# testing attenuation
#=
np2 = nextpow(2, 2 * length(tgrid));
fnpow2grid = FFTW.fftfreq(np2, inv(step(tgrid)));
fill!(medium.m[:Q], 100)
medium.fc = NamedArray([20.0, 100 * 2 * pi], ([:freqmin, :freqmax],))
medium.ic = NamedArray([5], ([:nsls],))

vp0 = medium.ref[:vp]
rho0 = medium.ref[:rho]
Ku = rho0 * vp0 * vp0
tau_epsilon = medium[:tau_epsilon]
tau_sigma = medium[:tau_sigma]

Kc = GeoPhyInv.complexK(Ku, 2.0 * pi .* fnpow2grid, tau_sigma, tau_epsilon)
=#






function check()
    # least-squares misfit
    paerr = GeoPhyInv.VNamedD_misfit(rec1, pa.c.data[1])
    err = GeoPhyInv.func_grad!(paerr)

    # normalize error
    error = err[1] / paerr.ynorm

    # desired accuracy?
    @test error < 1e-2
end


medium = Medium(:elastic_homo3D);
ageom=AGeom(medium.grid, SSrcs(2), [Srcs(3), Srcs(2)], [Recs(10), Recs(20)])
wav, tgrid=ricker(medium, 4, 1, 0.4)
srcwav = Srcs(tgrid, ageom, [:p])
update!(srcwav, [:p], wav)



pa = SeisForwExpt(
    FdtdElastic(),
    npw = 1,
    medium = medium,
    ageom = [ageom],
    srcwav = [srcwav],
    sflags = [2],
    rflags = [1],
    tgrid = tgrid,
    verbose = true,
);

@time update!(pa);




vp0 = mean(medium[:vp])
rho0 = mean(medium[:rho])
rec1 = GeoPhyInv.Born.mod(
    medium,
    medium_pert = medium,
    ageom = ageom,
    srcwav = srcwav,
    tgridmod = tgrid,
    src_flag = 2,
)



check()

