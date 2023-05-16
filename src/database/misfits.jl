
function lossvalue(loss, dataobs::NamedD{T}, data::NamedD{T}) where T
	return mapreduce(+, dataobs.d, data.d) do d1, d2 # loop over fields
		mapreduce(+, d1, d2) do dd1, dd2 # loop over time, receiver
			loss(dd1, dd2)
		end
  	end
end
function lossvalue(loss, dataobs::VNamedD, data::VNamedD)
	return mapreduce(+, dataobs, data) do d1, d2 # loop over supersources
		lossvalue(loss, d1, d2)
  	end
end


function gradient!(buffer, loss, dataobs::NamedD{T}, data::NamedD{T}) where T
	map(buffer.d, dataobs.d, data.d) do g, d1, d2 # loop over fields
		@. g = deriv(loss, d1, d2)
  	end
	return buffer
end
function gradient!(buffer, loss, dataobs::VNamedD, data::VNamedD)
	map(buffer, dataobs, data) do g, d1, d2 # loop over supersources
		gradient!(g, loss, d1, d2)
  	end
end

