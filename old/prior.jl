
"""
Update prior in `PFWI` either using mediumm or mediumi

* `priori` : prior model on inversion grid
"""
function update_prior!(pa; priori=nothing, priorm=nothing)

	if(priori===nothing) # use current mediumi as prior!
		if(priorm===nothing)
			copyto!(pa.priori, pa.mediumi)
		else
			interp_spray!(priorm, pa.priori, :interp)
			copyto!(pa.priori, pa.mediumi)
		end
	else
		!(isequal(priori.grid,pa.mediumi.grid)) && error("priori model has to be on igrid")
		copyto!(pa.priori, priori)
	end
	# update pa.mx.prior
	copyto!(pa.mx.prior,pa.priori,pa.parameterization)
end
