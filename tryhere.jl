using DifferentialEquations
using PyPlot
using PyCall

np = pyimport("numpy")
plt = pyimport("matplotlib.pyplot")
# using Plots;gr()

function f(dy,u,p,t)
    dy[1] = t.*t#(u[1])^2
    return dy
end

u0 = [0.1]
tspan = (0.1,4.0)
p= 2.0
prob = ODEProblem{true}(f,u0,tspan,p)
sol = solve(prob, TRBDF2(),abstol=1e-9,reltol=1e-15)

yu = ((1/3) * (sol.t).^3)
# plt.plot(sol.t,sol.u)
# plt.plot(sol.t,(1/3).*(sol.t).^3)
plt.plot(sol.t,np.array(sol.u) .- yu)

display(gcf())
plt.close()
