using Optim
initial_x = zeros(10)
buffer = rand(10) # Preallocate an appropriate buffer
last_x = similar(initial_x)
df = OnceDifferentiable(x -> f(x, buffer, initial_x),
			 (x, stor) -> g!(x, stor, buffer, last_x))
					 optimize(df, initial_x)



function calculate_common!(x, last_x, buffer)
if x != last_x
copy!(last_x, x)
#do whatever common calculations and save to buffer
println("common calculations")
end
end

function f(x, buffer, last_x)
calculate_common!(x, last_x, buffer)
x = buffer +1 # depends on buffer
return norm(x)
end

function g!(x, stor, buffer, last_x)
calculate_common!(x, last_x, buffer)
stor = buffer +1
end
