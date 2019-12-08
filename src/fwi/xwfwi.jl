"""
Update pa.w
"""
function wfwi!(pa::FWI; store_trace::Bool=true, extended_trace::Bool=false, time_limit=Float64=2.0*60., 
	       f_tol::Float64=1e-8, g_tol::Float64=1e-8, x_tol::Float64=1e-8)

	# convert initial model to the inversion variable
	x = zeros(wfwi_ninv(pa));
	last_x = rand(size(x)) # reset last_x

	println("updating w...")
	println("> wfwi: number of inversion variables:\t", length(x)) 

	# initial w to x
	Coupling_x!(x, pa, 1)

	(iszero(pa.paTD.y)) && error("dcal zero wfwi")

	f=x->func_grad_Coupling!(nothing, x,  pa)
	g! =(storage,x)->func_grad_Coupling!(storage, x,  pa)

	"""
	Unbounded LBFGS inversion, only for testing
	"""
	res = optimize(f, g!, x, 
		       LBFGS(),
		       Optim.Options(g_tol = g_tol, f_tol=f_tol, x_tol=x_tol,
		       iterations = 1000, store_trace = store_trace,
		       extended_trace=extended_trace, show_trace = true))
	"testing gradient using auto-differentiation"
#	res = optimize(f, x, 
#		       LBFGS(),
#		       Optim.Options(g_tol = 1e-12,
#	 	       iterations = 10, store_trace = true,
#		       extended_trace=true, show_trace = true))

	pa.verbose ? println(res) : nothing
	# update w in pa

	Coupling_x!(Optim.minimizer(res), pa, -1)

	f = Optim.minimum(res)
	return f
end # wfwi


#=

function xwfwi!(pa; max_roundtrips=100, max_reroundtrips=10, FWIAM_func=nothing, roundtrip_tol=1e-6,
		     optim_tols=[1e-6, 1e-6])

	if(FWIAM_func===nothing)
		FWIAM_func=x->FWIAM(x, optim_tols=optim_tols,name="FWI with Source Wavelet Estimation",
				    roundtrip_tol=roundtrip_tol, max_roundtrips=max_roundtrips,
				    max_reroundtrips=max_reroundtrips,
				    min_roundtrips=10,
				    reinit_func=x->initialize!(pa),
				    after_reroundtrip_func=x->(err!(pa); update_calsave!(pa);),
				    )
	end

	
	# create alternating minimization parameters
	f1=x->FWI.wfwi!(pa, extended_trace=false)
	f2=x->FWI.xfwi!(pa, extended_trace=false)
	paam=FWIAM_func([f1, f2])

	# do inversion
	go(paam)

	# print errors
	err!(pa)
	println(" ")
end

=#

