

"""
Returns the variance of data

"""
function Statistics.var(data1::TD)
	σ=0.0
	μ=Statistics.mean(data1)
	n=0
	for ifield = 1:length(data1.fields), iss = 1:data1.ageom.nss 
		for ir = 1:data1.ageom.nr[iss], it = 1:length(data1.tgrid)
			n += 1
			σ += (data1.d[iss, ifield][it, ir]-μ)^2 
		end
	end
	return σ*inv(n)
end

function Statistics.mean(data1::TD)
	n=0
	μ=0.0
	for ifield = 1:length(data1.fields), iss = 1:data1.ageom.nss 
		for ir = 1:data1.ageom.nr[iss], it = 1:length(data1.tgrid)
			n += 1
			μ += data1.d[iss, ifield][it, ir]
		end
	end
	return μ*inv(n)
end

