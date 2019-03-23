```@meta
EditURL = "https://github.com/TRAVIS_REPO_SLUG/blob/master/"
```

he linearized forward modeling operator `F`, such that
`Fx` can be computed without explicitly storing the operator matrix (
 see `LinearMaps.jl`);
the imaging/migration operator `F*`;

These maps are the building blocks of iterative optimization schemes.

```@example born_map
for scenario in [:downhole, :pizza]
	println("@@@@@@@@@@@@TESTING ", scenario)
	for rfields in [[:P], [:Vx], [:Vz]]
		pa=JG.SeisInvExpt(scenario, born_flag=true, rfields=rfields)


		F=JF.operator_Born(pa);

		x1=randn(size(F,2))
		x2=randn(size(F,2))
		x12=x1.+x2


		d12=F*x12
		d1=F*x1
		δmodtt1=copy(pa.paf.c.δmodtt)
		d2=F*x2
		δmodtt2=copy(pa.paf.c.δmodtt)


		d12new=d1.+d2

		f=Misfits.error_squared_euclidean!(nothing, d12, d12new, nothing, norm_flag=true)

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


		for i in 1:3
			adjtest()
		end
	end
end
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

