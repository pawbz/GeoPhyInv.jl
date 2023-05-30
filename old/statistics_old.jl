

"""
Returns the variance of data
"""
function Statistics.var(data1::VNamedD)
	σ=0.0
	μ=Statistics.mean(data1)
	n=0
	for iss in 1:length(data1)
		for ifield in names(data1[iss].d)[1]
			for ir = 1:data1[iss].n, it = 1:length(data1[iss].grid)
				n += 1
				σ += (data1[iss].d[ifield][it, ir]-μ)^2 
			end
		end
	end
	return σ*inv(n)
end

function Statistics.mean(data1::VNamedD)
	n=0
	μ=0.0
	for iss in 1:length(data1)
		for ifield in names(data1[iss].d)[1]
			for ir = 1:data1[iss].n, it = 1:length(data1[iss].grid)
				n += 1
				μ += data1[iss].d[ifield][it, ir]
			end
		end
	end
	return μ*inv(n)
end


function addnoise!(dataN::VNamedD, data::VNamedD, SNR)

	σx=Statistics.var(data)

	σxN=sqrt(σx^2*inv(10^(SNR/10.)))
	
	# factor to be multiplied to each scalar
	α=sqrt(σxN)
	for iss in 1:length(data)
		for ifield in names(data[iss].d)[1]
			for ir = 1:data[iss].n, it = 1:length(data[iss].grid)
				dataN[iss].d[ifield][it, ir] = dataN[iss].d[ifield][it, ir] + α*Random.randn()
			end
		end
	end
end



