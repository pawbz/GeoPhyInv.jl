
using Optim


prob = Optim.UnconstrainedProblems.examples["Rosenbrock"]


df = OnceDifferentiable(x -> prob.f(x), (x, storage) -> prob.g!(x, storage))

# without preconditioner
res_hz = optimize(df, prob.initial_x, [-1., -1.0], [2., 2.,], Fminbox();   optimizer=LBFGS, optimizer_o=Optim.Options(show_trace=true))
println(res_hz)


# testing preconditioner
P=spdiagm([1.,0.003], (0),2,2);
res_hz = optimize(df, prob.initial_x, [-1., -1.0], [2., 2.,], Fminbox(); P=P, optimizer=LBFGS, optimizer_o=Optim.Options(show_trace=true))
println(res_hz)

