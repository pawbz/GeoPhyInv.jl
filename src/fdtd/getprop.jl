"""
Return snaps stored in `SeisForwExpt` after using `mod!`. 

# Arguments

* `iss::Int64=1` : supersource index
"""
function Base.getindex(pa::PFdtd, s::Symbol, iss::Int=1)
	@assert iss≤length(pa.c.ageom[1])
	@assert s in [:snaps]
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
	end
end


