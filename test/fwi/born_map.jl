
using Test
using GeoPhyInv
using LinearMaps


for scenario in [:downhole, :pizza]
	println("@@@@@@@@@@@@TESTING ", scenario)
	for rfields in [[:p], [:vx], [:vz]]
		pa=SeisInvExpt(FdtdAcoustic{Born}(), LS(), scenario, rfields=rfields)


		F=LinearMap(pa);

		x1=randn(size(F,2)) 
		x2=randn(size(F,2)) 
		x12=x1.+x2


		d12=F*x12
		d1=F*x1
		δmodKI1=copy(pa.paf.c.δmod[:invK])
		d2=F*x2
		δmodKI2=copy(pa.paf.c.δmod[:invK])


		d12new=d1.+d2

		f=squared_euclidean!(nothing, d12, d12new, nothing, norm_flag=true)

		@test f<1e-25


		function adjtest()
			x=randn(size(F,2))
			y=randn(size(F,1))
			a=LinearAlgebra.dot(y,F*x)
			b=LinearAlgebra.dot(x,adjoint(F)*y)
			c=LinearAlgebra.dot(x, transpose(F)*F*x)
			println("adjoint test: ", a, "\t", b)       
			@test isapprox(a,b,rtol=1e-5)
			println("must be positive: ", c)
			@test c>0.0
		end


		for i in 1:1
			adjtest()
		end
	end
end




