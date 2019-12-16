using GeoPhyInv

medium_true=Medium(:pizza)

pa=SeisInvExpt(Fdtd(), LS(), :pizza);

update!(pa,solver=:ipopt);

