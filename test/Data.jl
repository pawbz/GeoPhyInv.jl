using Revise
using JuMIT
using Base.Test
using BenchmarkTools
using Calculus


fields=[:P]
# testing misfit and gradient calculation for TD objects
tgrid1=Grid.M1D(0.,1.,10);

tgrid2=Grid.M1D(0.,1.,20);

acqgeom=JuMIT.Acquisition.Geom_fixed(10,10,10,10,10,10,1,10)

x=JuMIT.Data.TD_ones(fields,tgrid1,acqgeom)
y=JuMIT.Data.TD_ones(fields,tgrid2,acqgeom)

# generate some random data
randn!(y.d[1,1])
randn!(x.d[1,1])

for func_attrib in [:cls, :xcorrcls]
	coup=JuMIT.Coupling.TD_delta(y.tgrid, [0.1,0.1], 0.0,  x.fields, x.acqgeom)
	randn!(coup.ssf[1,1])
	pa=JuMIT.Data.P_misfit(x,y, func_attrib=func_attrib, coup=coup);


	function err(xvec)

		for i in eachindex(x.d[1,1])
			pa.x.d[1,1][i]=xvec[i]
		end

		return JuMIT.Data.func_grad!(pa)
	end


	xvec=vec(pa.x.d[1,1])
	gg2=Calculus.gradient(x -> err(x), xvec)

	@time JuMIT.Data.func_grad!(pa,:dJx);
	gg1=vec(pa.dJx.d[1,1])
	# check gradient with Finite Differencing
	@test gg1 ≈ gg2

	function errw(xvec)

		for i in eachindex(pa.coup.ssf[1,1])
			pa.coup.ssf[1,1][i]=xvec[i]
		end

		return JuMIT.Data.func_grad!(pa)
	end

	xvec=vec(pa.coup.ssf[1,1])
	gg2=Calculus.gradient(x -> errw(x), xvec)

	@time JuMIT.Data.func_grad!(pa,:dJssf);
	gg1=vec(pa.dJssf[1,1])
	# check gradient with Finite Differencing
	@test gg1 ≈ gg2


end



