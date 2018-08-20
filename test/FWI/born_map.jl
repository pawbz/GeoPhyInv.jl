
pa, model=JG.xfwi_problem(:pizza, born_flag=true)


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

@test f<1e-30


function adjtest()
	x=randn(size(F,2))
	y=randn(size(F,1))
	a=LinearAlgebra.dot(y,F*x)
	b=LinearAlgebra.dot(x,adjoint(F)*y)
	c=LinearAlgebra.dot(x, transpose(F)*F*x)
	println("adjoint test: ", a, "\t", b)       
	@test isapprox(a,b,rtol=1e-6)
	println("must be positive: ", c)
	@test c>0.0
end


for i in 1:5
	adjtest()
end




