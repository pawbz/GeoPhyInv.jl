__precompile__()

module IO

import JuMIT.Data

type SegyHeader
	tracl::Int32
	tracr::Int32
	fldr::Int32
	tracf::Int32
	ep::Int32
	cdp::Int32
	cdpt::Int32
	trid::Int16
	nva::Int16
	nhs::Int16
	duse::Int16
	offset::Int32
	gelev::Int32
	selev::Int32
	sdepth::Int32
	gdel::Int32
	sdel::Int32
	swdep::Int32
	gwdep::Int32
	scalel::Int16
	scalco::Int16
	sx::Int32
	sy::Int32
	gx::Int32
	gy::Int32
	counit::Int16
	wevel::Int16
	swevel::Int16
	sut::Int16
	gut::Int16
	sstat::Int16
	gstat::Int16
	tstat::Int16
	laga::Int16
	lagb::Int16
	delrt::Int16
	muts::Int16
	mute::Int16
	ns::Int16
	dt::Int16
	gain::Int16
	igc::Int16
	igi::Int16
	corr::Int16
	sfs::Int16
	sfe::Int16
	slen::Int16
	styp::Int16
	stas::Int16
	stae::Int16
	tatyp::Int16
	afilf::Int16
	afils::Int16
	nofilf::Int16
	nofils::Int16
	lcf::Int16
	hcf::Int16
	lcs::Int16
	hcs::Int16
	year::Int16
	day::Int16
	hour::Int16
	minute::Int16
	sec::Int16
	timbas::Int16
	trwf::Int16
	grnors::Int16
	grnofr::Int16
	grnlof::Int16
	gaps::Int16
	otrav::Int16
	d1::Float32
	f1::Float32
	d2::Float32
	f2::Float32
	ungpow::Float32
	unscale::Float32
	ntr::Int32
	mark::Int16
	unass::Int16
end

segy_count = Dict{AbstractString,Int32}()
segy_count["tracl"]  = 0
segy_count["tracr"]  = 4
segy_count["fldr"]   = 8
segy_count["tracf"]  = 12
segy_count["ep"]     = 16
segy_count["cdp"]    = 20
segy_count["cdpt"]   = 24
segy_count["trid"]   = 28
segy_count["nva"]    = 30
segy_count["nhs"]    = 32
segy_count["duse"]   = 34
segy_count["offset"] = 36
segy_count["gelev"]  = 40
segy_count["selev"]  = 44
segy_count["sdepth"] = 48
segy_count["gdel"]   = 52
segy_count["sdel"]   = 56
segy_count["swdep"]  = 60
segy_count["gwdep"]  = 64
segy_count["scalel"] = 68
segy_count["scalco"] = 70
segy_count["sx"]     = 72
segy_count["sy"]     = 76
segy_count["gx"]     = 80
segy_count["gy"]     = 84
segy_count["counit"] = 88
segy_count["wevel"]  = 90
segy_count["swevel"] = 92
segy_count["sut"]    = 94
segy_count["gut"]    = 96
segy_count["sstat"]  = 98
segy_count["gstat"]  = 100
segy_count["tstat"]  = 102
segy_count["laga"]   = 104
segy_count["lagb"]   = 106
segy_count["delrt"]  = 108
segy_count["muts"]   = 110
segy_count["mute"]   = 112
segy_count["ns"]     = 114
segy_count["dt"]     = 116
segy_count["gain"]   = 118
segy_count["igc"]    = 120
segy_count["igi"]    = 122
segy_count["corr"]   = 124
segy_count["sfs"]    = 126
segy_count["sfe"]    = 128
segy_count["slen"]   = 130
segy_count["styp"]   = 132
segy_count["stas"]   = 134
segy_count["stae"]   = 136
segy_count["tatyp"]  = 138
segy_count["afilf"]  = 140
segy_count["afils"]  = 142
segy_count["nofilf"] = 144
segy_count["nofils"] = 146
segy_count["lcf"]    = 148
segy_count["hcf"]    = 150
segy_count["lcs"]    = 152
segy_count["hcs"]    = 154
segy_count["year"]   = 156
segy_count["day"]    = 158
segy_count["hour"]   = 160
segy_count["minute"] = 162
segy_count["sec"]    = 164
segy_count["timbas"] = 166
segy_count["trwf"]   = 168
segy_count["grnors"] = 170
segy_count["grnofr"] = 172
segy_count["grnlof"] = 174
segy_count["gaps"]   = 176
segy_count["otrav"]  = 178
segy_count["d1"]     = 180
segy_count["f1"]     = 184
segy_count["d2"]     = 188
segy_count["f2"]     = 192
segy_count["ungpow"] = 196
segy_count["unscale"]= 200
segy_count["ntr"]    = 204
segy_count["mark"]   = 208
segy_count["unass"]  = 210
segy_count["trace"]  = 240

function InitSegyHeader()
	h = SegyHeader(0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0)
	return h
end

function GrabSegyHeader(stream,swap_bytes,nt,file_header_size,j)
	position = file_header_size + (240+nt*4)*(j-1)
	seek(stream,position)
	h = InitSegyHeader()
	if (swap_bytes == false)
		for field in fieldnames(h)
			setfield!(h, field, read(stream,typeof(getfield(h,field))) )
		end
	else
		for field in fieldnames(h)
			setfield!(h, field, bswap(read(stream,typeof(getfield(h,field)))) )
		end
	end

	return h
end

function PutSegyHeader(stream,h,nt,file_header_size,j)
	position = file_header_size + (240+nt*4)*(j-1)
	seek(stream,position)
	for field in fieldnames(h)
		write(stream,getfield(h,field))
	end
end


"""
Read an input SeismicUnix file and output 
the time-domain data as `Data.TD`.

# Arguments
* `fname::AbstractString` : name of the SeismicUnix file
TODO: conversion to `TD`, currently the ouput is just a data matrix.
"""
function readsu_data(;
		     fname::AbstractString="", verbose::Bool=false
		     )
end


"""
Read a SU file and return trace headers along with data matrix
TODO: use bswap accordingly?
      extend it to SEGY format as well
"""
function readsu(fname, verbose=false)

	file_hsize=0
	swap_bytes = false

	stream = open(fname)
	seek(stream, segy_count["ns"]+file_hsize)
	nt = read(stream,Int16)
	total = 60 + nt
	nrecords = round(Int,(filesize(stream)-file_hsize)/4/total)
	verbose && println("number of traces: ",nrecords)
	verbose && println("number of samples per trace: ",nt)

	h_segy = Array{SegyHeader}(nrecords)
	seek(stream,file_hsize + segy_count["trace"])
	h_segy[1] = GrabSegyHeader(stream,swap_bytes,nt,file_hsize,1)
	dt = h_segy[1].dt/1000000 # convert to seconds


	d=zeros(nt, nrecords)
	# reading data and headers
	for j=1:nrecords
		position = file_hsize + total*(j-1)*4 + segy_count["trace"]
		# goto trace position
		seek(stream,position)
		# read trace
		d[:,j] = convert(Array{Float64,1}, read(stream,Float32,nt))
		# read header, seek is inside this function
		h_segy[j] = GrabSegyHeader(stream,swap_bytes,nt,file_hsize,j)
	end
	close(stream)

	return reshape(d, (nt, nrecords)), h_segy

end

"""
Write a SU file using data and headers

# Arguments

* `fname` : filename
* `d` : data matrix, output of readsu, for example
* `h_segy` : header vector, see output of readsu
"""
function writesu(fname, d, h_segy::Vector{SegyHeader}=[InitSegyHeader() for irec in 1:size(d,2)])

	nt = size(d, 1)
	total = 60 + nt
	#nrecords = round(Int,(filesize(stream)-file_hsize)/4/total)
	nrecords = size(d, 2)

	file_hsize=0
	swap_bytes = false

	stream = open(fname, "w")

	# reading data and headers
	for j=1:nrecords
		# add nt in header if absent
		(h_segy[j].ns == 0) && (h_segy[j].ns = nt)


		PutSegyHeader(stream,h_segy[j],nt,file_hsize,j)

		# write down the trace
		dvecsp = convert(Vector{Float32}, d[:,j])
		position = file_hsize + total*(j-1)*4 + segy_count["trace"]
		seek(stream,position)
		write(stream, dvecsp)
	end
	close(stream)
end

function TD(d::Data.TD, attrib=:su)
	
	nrecords = sum(d.acqgeom.nr)
	nss = d.acqgeom.nss

	# create a data matrix of first field
	dp = hcat(d.d[collect(1:nss),1][:,:]...);

	# headers
	h_segy = [InitSegyHeader() for irec in 1:nrecords]
	irec = 0
	for iss=1:nss
		for ir=1:d.acqgeom.nr[iss]
			irec += 1
			h_segy[irec].ns = d.tgrid.nx
			h_segy[irec].dt = d.tgrid.δx*1000000 
		end
	end
	writesu(fname, d, h_segy)
end

end # module
