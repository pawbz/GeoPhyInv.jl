```@meta
EditURL = "@__REPO_ROOT_URL__/"
```

```@example media
using BenchmarkTools
using GeoPhyInv
using Test
```

some mem tests

some model

```@example media
mgrid = [range(0.0, stop=10.,step=0.05), range(0.0, stop=10.,step=0.05)];
model = Medium(mgrid);


vpb=[2100.,2200.];vsb=[1500, 1700]; rhob=[2100., 2300.]
update!(model, [:vp,:vs,:rho], [vpb, vsb, rhob]);
fill!(model)


update!(model, [:vp,:rho], randn_perc=1e-3)

nznx=prod(length.(model.mgrid))
```

test copyto!

```@example media
model0=similar(model);
update!(model0, model)
fill!(model0)
update!(model0, [:vp,:rho], randn_perc=1e-3)

@btime copyto!(model0, model)
@test isequal(model0, model)
```

