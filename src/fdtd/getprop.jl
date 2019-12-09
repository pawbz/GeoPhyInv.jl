"""
Return snaps stored in `SeisForwExpt` after using `mod!`. 

# Arguments

* `iss::Int64=1` : supersource index
"""
function Base.getindex(pa::PFdtd, s::Symbol, iss::Int=1)
	@assert iss≤length(pa.c.ageom[1])
	@assert s in [:snaps,:data, :ageom, :srcwav]
	if(s==:snaps)
		nzd,nxd=length.(pa.c.model.mgrid)
		snaps_all=zeros(nzd,nxd,length(pa.c.itsnaps))
		ip=findall(in.(iss,pa.sschunks))[1]
		issp=findall(x->x==iss,pa.sschunks[ip])[]
		@sync begin
			p=procs(pa.p)[ip]
			@sync remotecall_wait(p) do 
				pap=localpart(pa.p)
				for i in eachindex(snaps_all)
					snaps_all[i]=pap.ss[issp].snaps[i]
				end
			end
		end
		return snaps_all
	elseif(s==:data)
		(iss==0) ?  (return pa.c.data[1]) :  (return pa.c.data[1][iss])
	elseif(s==:ageom)
		return pa.c.ageom
	elseif(s==:srcwav)
		return pa.c.srcwav
	end
end

function Base.show(io::Base.IO, pa::PFdtd)
	println(typeof(pa), "")
	println("pa[:snaps]")
	println("pa[:data] : modeled data after running `update!`")
	println("pa[:ageom]")
	println("pa[:srcwav]")
end

Base.show(io::Base.IO, pa::PFdtdp)=nothing
Base.show(io::Base.IO, pa::PFdtdc)=nothing
Base.show(io::Base.IO, pa::PFdtdss)=nothing



