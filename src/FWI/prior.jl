
"""
Update prior in `Param` either using modm or modi

* `priori` : prior model on inversion grid
"""
function update_prior!(pa; priori=nothing, priorm=nothing)

	if(priori===nothing) # use current modi as prior!
		if(priorm===nothing)
			copyto!(pa.priori, pa.modi)
		else
			Models.interp_spray!(priorm, pa.priori, :interp)
			copyto!(pa.priori, pa.modi)
		end
	else
		!(isapprox(priori,pa.modi)) && error("priori model has to be on igrid")
		copyto!(pa.priori, priori)
	end
	# update pa.mx.prior
	Models.Seismic_get!(pa.mx.prior,pa.priori,pa.parameterization)
end
