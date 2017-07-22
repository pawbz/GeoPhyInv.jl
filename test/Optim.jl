
using Optim



rosenbrock(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

df = OnceDifferentiable(x -> rosenbrock(x))

# without preconditioner
res_hz = optimize(df, [0., 0.], [-1., -1.0], [2., 2.,], Fminbox())
#res_hz = optimize(df, prob.initial_x) 
#		  optimizer_o=Optim.Options(show_trace=true))
println(res_hz)


# testing preconditioner
#P=spdiagm([1.,0.003], (0),2,2);
#res_hz = optimize(df, prob.initial_x, [-1., -1.0], [2., 2.,], Fminbox(); P=P, optimizer=LBFGS, optimizer_o=Optim.Options(show_trace=true))
#println(res_hz)

