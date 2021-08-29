using GeoPhyInv

medium_true=Medium(:pizza)

pa=SeisInvExpt(FdtdAcou(), LS(), :pizza);

update!(pa,solver=:ipopt);

