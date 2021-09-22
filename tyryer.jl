using DifferentialEquations
using Plots;gr()

jimmy = 0.5
function f(du,u,p,t)
    du[1,1] =p[1]*u[1,1]
    du[1,2] =p[2]*u[1,2]^1.2
    du[2,1] =0.5*p[3]*u[2,1]^1.3
    du[2,2] =p[4]*u[2,2]^1.4
end

#f(u,p,t) = 1.1*u
u0 = [1.0 1.0 
    2.0 2.0]
tspan = (0.1,1.0)
p=(0.2,0.3,0.4,0.5)
prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob, alg_hints=[:stiff],progress=true)#, reltol=1e-8, abstol=1e-8)
#sol = solve(prob, Tsit5(),dtmin=10^-17,maxiters = 5*10^8) #, reltol=1e-8, abstol=1e-8)

plot(sol)
#typeof(sol.t)