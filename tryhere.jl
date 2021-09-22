using DifferentialEquations
# using Plots;gr()

function f(dy,u,p,t)
    dy[1] = p*t*t#(u[1])
    dy[1] = dy[1]/10
    return dy
end

u0 = [0.0]
tspan = (0.0,4.0)
p= 2.0
prob = ODEProblem{true}(f,u0,tspan,p)
sol = solve(prob, alg_hints=[:stiff],progress=true)

plt.plot(sol.t,sol.u)
display(gcf())
plt.close()
