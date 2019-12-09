using GeoPhyInv
using Base.Test
using BenchmarkTools
using Calculus

const JD=GeoPhyInv.Data


fields=[:P]
# testing misfit and gradient calculation for TD objects
tgrid1=range(0., stop=1., length=5);
tgrid2=range(0., stop=1., length=10);

ageom=GeoPhyInv.Acquisition.AGeom_fixed(10,10,10,10,10,10,10,10)

x=GeoPhyInv.Data.TD_ones(fields,tgrid1,ageom)


yvec=randn(length(x))
xvec=similar(yvec)

@testset "testing copy!" begin
	@btime copy!(x, yvec);

	@btime copy!(xvec, x);

	@test xvec == yvec
end



@testset "simple LS error: x and y same time grid" begin
	y=GeoPhyInv.Data.TD_ones(fields,tgrid1,ageom)

	randn!(x)
	randn!(y)

	pa=GeoPhyInv.Data.P_misfit(x,y);
	@time GeoPhyInv.Data.func_grad!(pa,:dJx);
	gg1=vec(pa.dJx)


	xvec=vec(pa.x)
	gg2=Calculus.gradient(x -> (copy!(pa.x, x); return GeoPhyInv.Data.func_grad!(pa)), xvec)

	# check gradient with Finite Differencing
	@test gg1 ≈ gg2
end

#=


rrrrrr




# loop over same time grid and different time grid (interp_flag on/off)
for y in [GeoPhyInv.Data.TD_ones(fields,tgrid2,ageom), GeoPhyInv.Data.TD_ones(fields,tgrid1,ageom)]
	println("#########################################")


	# generate some random data
	randn!(y.d[1,1])
	randn!(x.d[1,1])

	for func_attrib in [:cls]
		coup=GeoPhyInv.Coupling.TD_delta(y.tgrid, [0.1,0.1], 0.0,  x.fields, x.ageom)
		randn!(coup.ssf[1,1])
		pa=GeoPhyInv.Data.P_misfit(x,y, func_attrib=func_attrib, coup=coup);



		xvec=vec(pa.x.d[1,1])
		gg2=Calculus.gradient(x -> err(x), xvec)

		@time GeoPhyInv.Data.func_grad!(pa,:dJx);
		gg1=vec(pa.dJx.d[1,1])
		# check gradient with Finite Differencing
		@test gg1 ≈ gg2

		function errw(xvec)

			for i in eachindex(pa.coup.ssf[1,1])
				pa.coup.ssf[1,1][i]=xvec[i]
			end

			return GeoPhyInv.Data.func_grad!(pa)
		end

		xvec=vec(pa.coup.ssf[1,1])
		gg2=Calculus.gradient(x -> errw(x), xvec)

		@time GeoPhyInv.Data.func_grad!(pa,:dJssf);
		gg1=vec(pa.dJssf[1,1])
		# check gradient with Finite Differencing
		@test gg1 ≈ gg2


	end
end



=#
