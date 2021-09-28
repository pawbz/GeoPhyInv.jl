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





# without attenuation 
medium = Medium(:acou_homo1, 5);
ageom = AGeom(medium.mgrid, :xwell);
tgrid = range(0.0, stop = 2.0, length = 2000)
wav = ricker(10.0, tgrid, tpeak = 0.25);
srcwav = SrcWav(tgrid, ageom, [:p])
update!(srcwav, [:p], wav)


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



pa = SeisForwExpt(
    FdtdAcou(),
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

check()


# @info "testing with attenuation" 
# medium = Medium(:acou_homo1);
# medium = Medium(medium.mgrid, [:vp, :rho, :Q])
# update!(medium, [:vp, :rho, :Q], [[1500, 2500], [1500, 2500], [10, 10]])
# fill!(medium)

# ageom = AGeom(medium.mgrid, :xwell);
# tgrid = range(0.0, stop = 2.0, length = 1000)
# wav = ricker(10.0, tgrid, tpeak = 0.25);
# srcwav = SrcWav(tgrid, ageom, [:p])
# update!(srcwav, [:p], wav)

# pa = SeisForwExpt(
#     FdtdAcouVisco(),
#     npw = 1,
#     medium = medium,
#     ageom = [ageom],
#     srcwav = [srcwav],
#     sflags = [2],
#     rflags = [1],
#     tgrid = tgrid,
#     verbose = true,
# );


# @time update!(pa);



# vp0 = mean(medium[:vp])
# rho0 = mean(medium[:rho])
# rec1 = GeoPhyInv.Born.mod(
#     pa.c.exmedium,
#     medium_pert = pa.c.exmedium,
#     ageom = ageom,
#     srcwav = srcwav,
#     tgridmod = tgrid,
#     src_flag = 2,
# )


# check()
