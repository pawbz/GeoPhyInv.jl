
pa, mod=JG.xfwi_problem(:pizza, born_flag=true)


F=JF.operator_Born(pa);

fac=1e-1
x1=randn(size(F,2)) .* fac
#fill!(x1,0.0)
#x1[10]=fac
x2=randn(size(F,2)) .* fac
#fill!(x2,0.0)
#x2[220]=-1.*fac
x12=x1.+x2


d12=F*x12
d1=F*x1
δmodtt1=copy(pa.paf.c.δmodtt)
d2=F*x2
δmodtt2=copy(pa.paf.c.δmodtt)


d12new=d1.+d2


f=Misfits.error_squared_euclidean!(nothing, d12, d12new, nothing, norm_flag=true)
println(d12-d1-d2)

println(f)






