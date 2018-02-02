using Revise
using JuMIT
using Base.Test


fields=[:P]
# testing misfit and gradient calculation for TD objects
tgrid1=JuMIT.Grid.M1D(0.,1.,10);

tgrid2=JuMIT.Grid.M1D(0.,1.,20);

acqgeom=JuMIT.Acquisition.Geom_fixed(10,10,10,10,10,10,1,10)

x=JuMIT.Data.TD_ones(fields,tgrid1,acqgeom)
y=JuMIT.Data.TD_ones(fields,tgrid2,acqgeom)
randn!(y.d[1,1])
randn!(x.d[1,1])

for func_attrib in [:cls, :xcorrcls]
	pa=JuMIT.Data.Param_error(x,y, func_attrib=func_attrib);


	function err(xvec)

		for i in eachindex(x.d[1,1])
			pa.x.d[1,1][i]=xvec[i]
		end

		return JuMIT.Data.error!(pa)
	end


	xvec=vec(pa.x.d[1,1])
	g2=zeros(xvec)
	JuMIT.Inversion.finite_difference!(x -> err(x), xvec, g2, :central)

	JuMIT.Data.error!(pa,:dJx);
	g1=vec(pa.dJx.d[1,1])
	@test g1 ≈ g2
end



