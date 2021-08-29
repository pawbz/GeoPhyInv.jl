using GeoPhyInv
#md using Plots; gr();

# ### Medium
medium_true=Medium(:pizza)


# ### SeisInvExpt
pa=SeisInvExpt(FdtdAcou(), LS(), :pizza);

# ### Initial plotting
#md p1=heatmap(pa.mediumm, :vp, title="Initial vp")
#md scatter!(pa.ageom, SSrcs())
#md scatter!(pa.ageom, Recs())
#md p2=heatmap(medium_true, :vp, title="True vp")
#md scatter!(pa.ageom, SSrcs())
#md scatter!(pa.ageom, Recs())
#md plot(p2, p1, size=(700,300))


# ### Final plotting 
update!(pa,solver=:ipopt);

# ### Plotting #1
#md p1=heatmap(pa.mediumm, :vp, title="Inverted vp")
#md scatter!(pa.ageom, SSrcs())
#md scatter!(pa.ageom, Recs())
#md p2=heatmap(medium_true, :vp, title="True vp")
#md scatter!(pa.ageom, SSrcs())
#md scatter!(pa.ageom, Recs())
#md plot(p2, p1, size=(700,300))




