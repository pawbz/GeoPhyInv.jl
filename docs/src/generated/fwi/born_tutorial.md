```@meta
EditURL = "<unknown>/fwi/born_tutorial.jl"
```

```@example born_tutorial
using Test
using GeoPhyInv
using LinearMaps
using LinearAlgebra
```

### SeisInvExpt

```@example born_tutorial
pa=SeisInvExpt(FdtdAcouBorn(), LS(), :pizza, rfields=[:P, :vz]); # use linearized modeling attribute
F=LinearMap(pa); # generate LinearMap
nothing #hide
```

### Linearity test

```@example born_tutorial
x1=randn(size(F,2)); # random input 1
x2=randn(size(F,2)); # random input 2
x12=x1.+x2; # sum
nothing #hide
```

```@example born_tutorial
d12=F*x12;
d1=F*x1; #
d2=F*x2;
d12new=d1.+d2;
nothing #hide
```

```@example born_tutorial
```

### Adjoint test

```@example born_tutorial
x=randn(size(F,2)); # random x
y=randn(size(F,1)); # random y
a=LinearAlgebra.dot(y,F*x);
b=LinearAlgebra.dot(x,adjoint(F)*y);
c=LinearAlgebra.dot(x, transpose(F)*F*x);
nothing #hide
```

```@example born_tutorial
println("adjoint test: ", a, "\t", b)
@test isapprox(a,b,rtol=1e-5)
```

```@example born_tutorial
println("x*transpose(F)*F*x be positive: ", c)
@test c>0.0
```

