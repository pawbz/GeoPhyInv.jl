```@meta
EditURL = "<unknown>/fwi/pizza.jl"
```

```@example pizza
using GeoPhyInv
using Plots; gr();
nothing #hide
```

### Medium

```@example pizza
medium_true=Medium(:pizza)
```

### SeisInvExpt

```@example pizza
pa=SeisInvExpt(FdtdAcou(), LS(), :pizza);
nothing #hide
```

### Initial plotting

```@example pizza
p1=heatmap(pa.mediumm, :vp, title="Initial vp")
scatter!(pa.ageom, SSrcs())
scatter!(pa.ageom, Recs())
p2=heatmap(medium_true, :vp, title="True vp")
scatter!(pa.ageom, SSrcs())
scatter!(pa.ageom, Recs())
plot(p2, p1, size=(700,300))
```

### Final plotting

```@example pizza
update!(pa,solver=:ipopt);
nothing #hide
```

### Plotting #1

```@example pizza
p1=heatmap(pa.mediumm, :vp, title="Inverted vp")
scatter!(pa.ageom, SSrcs())
scatter!(pa.ageom, Recs())
p2=heatmap(medium_true, :vp, title="True vp")
scatter!(pa.ageom, SSrcs())
scatter!(pa.ageom, Recs())
plot(p2, p1, size=(700,300))
```

