module IO

import SIT.F90libs


function readsu_data(;
		     fname::AbstractString=""
		     )

	nt = [0]; nrecords = [0];
	# read number of time samples and number of records first 
	ccall( (:readsu_nt_nrecords, F90libs.io), Void, 
       (Ptr{UInt8}, Ref{Int64}, Ref{Int64}),
	      fname, nt, nrecords); 
	nt = nt[1]; nrecords = nrecords[1];

	println("readsu_data: nt = ", string(nt))
	println("readsu_data: nrecords = ", string(nrecords))
	data = zeros(nt, nrecords);

	# read data now 
	ccall( (:readsu_data_fast, F90libs.io), Void, 
	      (Ptr{UInt8}, Ref{Int64}, Ref{Int64}, Ptr{Float64}),
	      fname, nt, nrecords, data);

	max(abs(data) == 0.0) && warn("readsu_data is zeros")
	return reshape(data, (nt, nrecords)), nt, nrecords
end



end # module
