using SIT
using Base.Test

model_pert = SIT.Gallery.Seismic(:acou_homo2);
SIT.Models.Seismic_addon!(model_pert, circ_rad=30., circ_loc=[0.,0.],circ_pert=0.2, fields=[:χvp])
model0 = SIT.Gallery.Seismic(:acou_homo2);

tgrid = SIT.Gallery.M1D(:acou_homo2);
acqgeom = SIT.Gallery.Geom(model0.mgrid,:oneonev);
acqsrc = SIT.Gallery.Src(:acou_homo2);


# Store boundaries
# Generate Born Data
buffer = SIT.Fdtd.mod(npropwav=2, model=model0, model_pert=model_pert, 
	acqgeom=[acqgeom,acqgeom], acqsrc=[acqsrc,acqsrc], 
	src_flags=[2.0, 0.0], recv_flags = [0.0, 2.0], 
	tgridmod=tgrid, verbose=true, boundary_save_flag=true, born_flag=true);

#rec1 = SIT.Fdtd.fdtd_born_mod();


	# source sinks
		acqsrcsink = deepcopy(acqsrc); 
		acqsrcsink.wav = [-1.0.*flipdim(acqsrc.wav[i,j],1) for 
			   i in 1:acqsrc.nss, j in 1:acqsrc.nfield];

		
# adjoint simulation
	adjsrc = SIT.Inversion.AdjSrc(buffer[1][1])
	adjacq = SIT.Inversion.AdjGeom(acqgeom);
# migrate Born data
		adj = SIT.Fdtd.mod(npropwav=2, model=model0,  
		     acqgeom=[acqgeom,adjacq], acqsrc=[acqsrcsink, adjsrc], src_flags=[2.0, -2.0], 
		     tgridmod=tgrid, grad_out_flag=true, boundary_in=buffer[2], verbose=true)


		g1=zeros(model0.mgrid.nz*model0.mgrid.nx)
		g2=zeros(model0.mgrid.nz*model0.mgrid.nx)
		SIT.Models.Seismic_chainrule!(adj[3],model_pert,g1,g2,[:χKI, :χρI],-1)


m1 = vec(SIT.Models.Seismic_get(model_pert, :χKI)-SIT.Models.Seismic_get(model0, :χKI))
m2 =  vec(SIT.Models.Seismic_get(model_pert, :χρI)-SIT.Models.Seismic_get(model0, :χρI))

println(maximum(m2))

gvec=vcat(g1, g2)
mvec=vcat(m1, m2)

# dot product test
@test_approx_eq SIT.Data.TD_dot(buffer[1][1], buffer[1][1]) dot(gvec, mvec) 


