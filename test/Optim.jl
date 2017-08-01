
using Optim, LineSearches


f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

function g!(storage, x)
	storage[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
	storage[2] = 200.0 * (x[2] - x[1]^2)
end

		#df = OnceDifferentiable(x -> f(x), (storage, x) -> g!(storage,x))
		df = OnceDifferentiable(f, g!)


# without preconditioner
res_hz = optimize(df, [0., 0.], [-1., -1.0], [2., 2.,], Fminbox{ConjugateGradient}(), 
		  linesearch = BackTracking(order=3), f_tol=1e-5, x_tol=1e-3,
		  iterations=100,
		  optimizer_o=Optim.Options( f_tol=1e-5, x_tol=1e-3,iterations=2, show_trace=true      ))
println(res_hz)

#res_hz = optimize(df, [0., 0.], LBFGS())
#res_hz = optimize(df, prob.initial_x) 
#		  optimizer_o=Optim.Options(show_trace=true))
#println(res_hz)


# testing preconditioner
#P=spdiagm([1.,0.003], (0),2,2);
#res_hz = optimize(df, prob.initial_x, [-1., -1.0], [2., 2.,], Fminbox(); P=P, optimizer=LBFGS, optimizer_o=Optim.Options(show_trace=true))
#println(res_hz)

