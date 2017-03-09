using SeismicInversion
using Plots

a, b, c = SeismicInversion.IO.readsu_data(fname="/home/pawbz/marmousi2/vp_marmousi-ii.su")
b, c = SeismicInversion.IO.readsu_data(fname="/home/pawbz/marmousi2/vp_marmousi-ii.su")
c = SeismicInversion.IO.readsu_data(fname="/home/pawbz/marmousi2/vp_marmousi-ii.su")

plot((SeismicInversion.Wavelets.ricker()))
plot((SeismicInversion.Wavelets.ricker(trim_tol=1e-6)))

println("acoustic modeling in homogeneous media...")

mesh = SeismicInversion.mesh.Mesh2D()

model = SeismicInversion.Models.Seismic("test_homo_acoustic")
model.vp0
model.vs0


acqgeom = SeismicInversion.Acquisition.Geom(2, 10);

acqgeom.nr
acqgeom.sz
records = SeismicInversion.Fdtd.fdtd_mod(acqgeom=acqgeom)
plot((records[:,1,1]))
heatmap(model.χvp)
gui()
mean(model.χvs)
heatmap(model.χvs)
close("all")

SeismicInversion.fdtd.fdtd_mod()
