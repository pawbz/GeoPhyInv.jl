
"""
Update prior in `Param` either using modm or modi

* `priori` : prior model on inversion grid
* `priorw` : prior weights on the inversion grid
"""
function update_prior!(pa; priori=nothing, priorw=nothing, priorm=nothing)

	if(priori===nothing) # use current modi as prior!
		if(priorm===nothing)
			copy!(pa.priori, pa.modi)
		else
			Models.interp_spray!(priorm, pa.priori, :interp)
			copy!(pa.priori, pa.modi)
		end
	else
		!(isapprox(priori,pa.modi)) && error("priori model has to be on igrid")
		copy!(pa.priori, priori)
	end
	# update pa.mx.prior
	Models.Seismic_get!(pa.mx.prior,pa.priori,pa.parameterization)


	if(!(priorw===nothing))
		!(isapprox(priorw,modi)) && error("priorw model has to be on igrid")
		copy!(pa.priorw, priorw)
		# update pa.mx.w
		Models.Seismic_get!(pa.mx.w,pa.priorw,pa.parameterization)
	end

end
