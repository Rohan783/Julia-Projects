using PyCall
using SpecialFunctions
using DifferentialEquations
using PyPlot
# using IRKGaussLegendre
# using LSODA
#using Plots;gr()

#scriptdir = @__DIR__
#pushfirst!(PyVector(pyimport("sys")."path"), scriptdir);
#gellmann = pyimport("gellmann")
np = pyimport("numpy")
sp = pyimport("sympy")
IPy = pyimport("IPython")
scipy = pyimport("scipy")
plt = pyimport("matplotlib.pyplot")
#Dagger = sp.physics.quantum.dagger.Dagger

#println(f.struc_cons())
#r=np.array(f.struc_cons())
#println(r)

delmass = np.array([0.1,1,2,4])
theta_all = np.array(asin.(sqrt.([3*10^-3,3*10^-2]))/2)

theta12 = deg2rad(33.45)
name="November17_TRBDF"
delta_sqm12 =  7.6*(10^(-5))
delta_sqm23 = delmass[1]
theta13 = theta_all[1]#deg2rad(0.642)
theta23 = theta_all[1]#deg2rad(0.5)#theta_all[1] 
del_cp = 0#np.radians(234)
# x = np.linspace(0.1,20,30)
Gf= 1.16*(10^(-23))

name = name*"_"*string(1+np.where(delmass.==delta_sqm23)[1][1])*string(1+np.where(theta_all.==theta13)[1][1])*string(1+np.where(theta_all.==theta23)[1][1])#*"_r15c"
M2_z= 91^2 #GeV #check code for powers of 10
G= 6.707*(10^-57)
cos2th_w= 0.7777
sec2th_w= 1/cos2th_w
#T = np.linspace(60*10^6,1*10^6,10000)

#f= gellmann.struc_cons()
#r=np.stack([f[1:end,1:end,1:end],f[1:end,1:end,1:end]])

# z=np.linspace(0,1,100)

# def fxy(y):
#     x_ext= 3.15
#     x_min= 10**-1
#     x_max= 50
#     return (y*x_ext*(x_max-x_min) + x_ext*x_min + x_max*x_min) /((1-y)*x_max + y*x_ext + x_min)
x= [3.15]
# #x=fxy(z)

W_cutoff = 4*sec2th_w

function Hubble(T)
    return (5.4*(T^2))/(1.22*(10^16))
end

function fermi(x,chem_p)
    return 1 / (1 .+ exp.(x.- chem_p))
end

f0 = fermi.(x,0)
f0x2 = f0.*(x.^2)
#norm = np.trapz(f0x2,x)
# # In[ ]:


m1 ,m2 ,m3 = sp.symbols("m1 m2 m3", real = true)
#th12,th23, th13 ,del_cp, delta_sqm12, delta_sqm23 ,p  = symbols('th12 th23 th13 del_cp d_12 d_23 p', real = True)
V= sp.ones(1,9)*0
s12= np.sin(theta12)
s13= np.sin(theta13)
s23= np.sin(theta23)
c12= np.cos(theta12)
c13= np.cos(theta13)
c23= np.cos(theta23)

U= sp.Matrix([[c12*c13, s12*c13, s13*sp.exp(-del_cp*1im)], 
           [-s12*c23 - c12*s23*s13*sp.exp(del_cp*1im), c12*c23 - s12*s23*s13*sp.exp(del_cp*1im), s23*c13], 
           [s12*s23 - c12*c23*s13*sp.exp(del_cp*1im), -c12*s23- s12*c23*s13*sp.exp(del_cp*1im), c23*c13]])
#display(U)

H_o= sp.diag(m1^2,m2^2,m3^2)
#display(H_o)

Ud= U.H
#display(Ud)
Hf = U*H_o*Ud
# display(Hf)

#f= gellmann.struc_cons()

V[1]=sp.simplify((Hf[1,1]+Hf[2,2]+Hf[3,3])/3)
V[2]=sp.re(Hf[1,2])
V[3]=sp.im(Hf[1,2])
V[4]=sp.collect((Hf[1,1]-Hf[2,2])/2,(m1 ,m2 ,m3))
V[5]=sp.re(Hf[1,3])
V[6]=sp.im(Hf[1,3])
V[7]=sp.re(Hf[2,3])
V[8]=sp.im(Hf[2,3])
V[9]=sp.collect(-np.sqrt(3)*(Hf[3,3]-V[1])/2,(m1 ,m2 ,m3))

function forward_v(T,Nu,aNu,L)
    a = (-1.5207*(10^-15)*(T^5)*x.*(Nu.+aNu))
    b = (3.996*(10^-6)*(T^3)*L)
    return [a .+ b,a .- b]
end

function forward_v4e(T,Nu,aNu,L)
    g_e = (1.0 .+ (W_cutoff)./(Nu .+ aNu))
    a = (-1.5207*(10^-15)*(T^5)*x.*(Nu.+aNu).*g_e)
    b = (3.996*(10^-6)*(T^3)*L)
    return [a .+ b,a .- b]
end

function cross12(T::Z,pol_diff1::Vector{Z},pol_diff2::Vector{Z}) where {Z <: Number}
    ans = 1.49851*(10^-6)*(T^3).*[pol_diff1,pol_diff2]
    return Matrix{(2,lenx),eltype(ans)}ans
end

sq3 = sqrt(3)
function GMcross(dpol_dt,pot_all,pol)
    # R_k = V_i*P_j*F_ijk
    # dpol_dtu=zeros(Real,18,lenx)
    #add +1 
    dpol_dt[3,:] = -pol[3+1,:].*pot_all[4-1,:] + pol[4+1,:].*pot_all[3-1,:] + 0.5*pol[6+1,:].*pot_all[7-1,:] + 0.5*pol[8+1,:].*pot_all[5-1,:]
    dpol_dt[4,:] = pol[2+1,:].*pot_all[4-1,:] - pol[4+1,:].*pot_all[2-1,:] - 0.5*pol[5+1,:].*pot_all[7-1,:] +0.5*pol[7+1,:].*pot_all[5-1,:]
    dpol_dt[5,:] = -pol[2+1,:].*pot_all[3-1,:] + pol[3+1,:].*pot_all[2-1,:] +0.5*pol[6+1,:].*pot_all[5-1,:] -0.5*pol[8+1,:].*pot_all[7-1,:]
    dpol_dt[6,:] = 0.5*pol[3+1,:].*pot_all[7-1,:] -0.5*pol[7+1,:].*pot_all[3-1,:] -0.5*pol[8+1,:].*pot_all[2-1,:] -0.5*(pot_all[4-1,:] + sq3*pot_all[9-1,:]).*pol[6+1,:]
    dpol_dt[7,:] = -0.5*pol[2+1,:].*pot_all[7-1,:] +0.5*pol[7+1,:].*pot_all[2-1,:] -0.5*pol[8+1,:].*pot_all[3-1,:] -0.5*(pol[4+1,:] .+ sq3*pol[9+1,:]).*pot_all[5-1,:] +0.5*(pot_all[4-1,:] .+ sq3*pot_all[9-1,:]).*pol[5+1,:]
    dpol_dt[8,:] = -0.5*pol[3+1,:].*pot_all[5-1,:] +0.5*pol[5+1,:].*pot_all[3-1,:] -0.5*pol[6+1,:].*pot_all[2-1,:] -0.5*(-pot_all[4-1,:] .+ sq3*pot_all[9-1,:]).*pol[8+1,:]
    dpol_dt[9,:] = -0.5*pol[2+1,:].*pot_all[5-1,:] +0.5*pol[5+1,:].*pot_all[2-1,:] +0.5*pol[6+1,:].*pot_all[3-1,:] -0.5*(-pol[4+1,:] .+sq3*pol[9+1,:]).*pot_all[7-1,:] +0.5*(-pot_all[4-1,:] .+ sq3*pot_all[9-1,:]).*pol[7+1,:]
    dpol_dt[10,:] = sq3*0.5*(pol[6+1,:] .* pot_all[5-1,:] .+ (pol[8+1,:] .* pot_all[7-1,:]))
    
    dpol_dt[3+8,:] = -(-pol[3+9,:].*pot_all[4+7,:] + pol[4+9,:].*pot_all[3+7,:] + 0.5*pol[6+9,:].*pot_all[7+7,:] +0.5*pol[8+9,:].*pot_all[5+7,:])
    dpol_dt[4+8,:] = -(pol[2+9,:].*pot_all[4+7,:] - pol[4+9,:].*pot_all[2+7,:] - 0.5*pol[5+9,:].*pot_all[7+7,:] +0.5*pol[7+9,:].*pot_all[5+7,:])
    dpol_dt[5+8,:] = -(-pol[2+9,:].*pot_all[3+7,:] + pol[3+9,:].*pot_all[2+7,:] +0.5*pol[6+9,:].*pot_all[5+7,:] -0.5*pol[8+9,:].*pot_all[7+7,:])
    dpol_dt[6+8,:] = -(0.5*pol[3+9,:].*pot_all[7+7,:] -0.5*pol[7+9,:].*pot_all[3+7,:] -0.5*pol[8+9,:].*pot_all[2+7,:] -0.5*(pot_all[4+7,:] +sq3*pot_all[9+7,:]).*pol[6+9,:])
    dpol_dt[7+8,:] = -(-0.5*pol[2+9,:].*pot_all[7+7,:] +0.5*pol[7+9,:].*pot_all[2+7,:] -0.5*pol[8+9,:].*pot_all[3+7,:] -0.5*(pol[4+9,:] +sq3*pol[9+9,:]).*pot_all[5+7,:] +0.5*(pot_all[4+7,:] +sq3*pot_all[9+7,:]).*pol[5+9,:])
    dpol_dt[8+8,:] = -(-0.5*pol[3+9,:].*pot_all[5+7,:] +0.5*pol[5+9,:].*pot_all[3+7,:] -0.5*pol[6+9,:].*pot_all[2+7,:] -0.5*(-pot_all[4+7,:] +sq3*pot_all[9+7,:]).*pol[8+9,:])
    dpol_dt[9+8,:] = -(-0.5*pol[2+9,:].*pot_all[5+7,:] +0.5*pol[5+9,:].*pot_all[2+7,:] +0.5*pol[6+9,:].*pot_all[3+7,:] -0.5*(-pol[4+9,:] +sq3*pol[9+9,:]).*pot_all[7+7,:] +0.5*(-pot_all[4+7,:] +sq3*pot_all[9+7,:]).*pol[7+9,:])
    dpol_dt[10+8,:] = -((sq3*0.5*pol[6+9,:].*pot_all[5+7,:]) .+ (sq3*0.5*pol[8+9,:].*pot_all[7+7,:]))
    # return dpol_dt[3:18,:]
end

coeff_del23 = convert(Vector{Float64},[sp.re(sp.simplify(V[i].coeff(m3^2))) for i=2:9])
coeff_del12 = convert(Vector{Float64},[sp.re(sp.simplify(V[i].coeff(m1^2))) for i=2:9])

vac_coeff = coeff_del23.*delta_sqm23 .- coeff_del12.*delta_sqm12
vac_for_all_x = repeat(Array(transpose(vac_coeff[[CartesianIndex()],:] ./ (2*x))),2)#np.einsum("i,j",vac_coeff,(0.5/x))

function vacuum_potential(T)
    return vac_for_all_x/(T*10^6)
end

# norm = 8*zeta(3)#np.trapz(fermi(x)*(x**3), x)
lenx=length(x)
GGG= (4.26984*10^-17) .*x
coeffs= repeat([0.5*4,0.5*2.9,1.97,0.5*4,0.5*2.9,1.97],1,2)
C_dash = repeat([1.23,1.23],1,2)
C_tau =0.92 ; C_mu = 0.92; C_e = 1.27; G_e = -0.3096 ; G_nu_tau = 0.51
v=zeros(Real,16,lenx)
#dpol_dt = zeros(Real,20,lenx)

twoxsqrt3 = 2*sqrt(3)
sqrt3 = sqrt(3)
root3by2 = sqrt(3)/2
zeta3= zeta(3)

Lscale = 1e-10
L_0 = 1e-10

function equation!(dpol_dt,pol::Matrix{Z},v,T) where {Z <: Number}
    T= (T/10^6)

    Nu_e= 0.5*(pol[1,:] .+ pol[5,:] .+ (pol[10,:]/sqrt3))
    Nu_mu= 0.5*(pol[1,:] .- pol[5,:] .+ (pol[10,:]/sqrt3))

    aNu_e= 0.5*(pol[2,:] .+ pol[13,:] .+ (pol[18,:]/sqrt3))
    aNu_mu= 0.5*(pol[2,:] .- pol[13,:] .+ (pol[18,:]/sqrt3))

    #Nu_e_norm, aNu_e_norm, Nu_mu_norm, aNu_mu_norm = ([Nu_e .* f0x2,aNu_e  .* f0x2,Nu_mu  .* f0x2,aNu_mu  .* f0x2])
    # L_e = 1e-10 .+ 2*(Nu_e - aNu_e) .+ (Nu_mu - aNu_mu)
    # L_mu = 1e-10 .+ 2*(Nu_mu - aNu_mu) .+ (Nu_e - aNu_e)
    # L_e = pol[19,:] .* Lscale
    # L_mu = pol[20,:] .* Lscale
    L_e =  L_0 .+ (3/8)*(2*pol[19,:] .+ pol[20,:]).*Lscale
    L_mu = L_0 .+ (3/8)*(2*pol[20,:] .+ pol[19,:]).*Lscale
    # chem_p_e = -3.628*pi*sinh.(asinh.(-1.2087*(L_e))/3)
    # chem_p_mu= -3.628*pi*sinh.(asinh.(-1.2087*(L_mu))/3)

    v[:,:] = vacuum_potential(T)

    f_veL,f_ve_L = forward_v4e(T,Nu_e,aNu_e,L_e)
    f_vmuL,f_vmu_L = forward_v(T,Nu_mu,aNu_mu,L_mu)

    cross_real, cross_imag = (1.49851*(10^-6)*(T^3)) .* [pol[3,:] .- pol[11,:], pol[4,:] .- pol[12,:]]

    #println(f_vmuL[3],f_vmu_L[3] ,cross_real[3], cross_imag[3])
    v[1,:] = v[1,:] .+ cross_real
    v[2,:] = v[2,:] .- cross_imag
    v[3,:]= v[3,:] .+ (f_veL - f_vmuL)/2
    v[8,:]= v[8,:] .+ (f_veL + f_vmuL)/twoxsqrt3

    v[9,:] = v[9,:] .- cross_real
    v[10,:]= -(v[10,:] .+ cross_imag)
    v[11,:]= v[11,:] .+ (f_ve_L - f_vmu_L)/2
    v[16,:]= v[16,:] .+ (f_ve_L + f_vmu_L)/twoxsqrt3

    GMcross(dpol_dt,2*v,pol) 
    # dpol_dt[19,:] .= 0
    # dpol_dt[20,:] .= 0
    dpol_dt[19,:] = ((dpol_dt[5,:]-dpol_dt[13,:]) + ((dpol_dt[10,:]-dpol_dt[18,:])*sqrt3))/(2*Lscale)
    dpol_dt[20,:] = ((dpol_dt[13,:]-dpol_dt[5,:]) + ((dpol_dt[10,:]-dpol_dt[18,:])*sqrt3))/(2*Lscale)

    G2T5 = - (GGG/3.15) * T^5

    R_e = 0.5*G2T5.*((Nu_e .-1) + (aNu_e .-1))
    R_mu = 0.5*G2T5.*(0.56*(Nu_mu .-1) + 0.56*(aNu_mu .-1))
    aR_e = 0.5*G2T5.*((aNu_e .-1) + (Nu_e .-1))
    aR_mu = 0.5*G2T5.*(0.56*(aNu_mu .-1) + 0.56*(Nu_mu .-1))

    dpol_dt[1,:] = (2/3)*(R_e .+ R_mu)
    dpol_dt[5,:] = dpol_dt[5,:] .+ (R_e .- R_mu)
    dpol_dt[10,:] = dpol_dt[10,:] .+ (dpol_dt[1,:] *root3by2)

    dpol_dt[2,:] = (2/3)*(aR_e .+ aR_mu)
    dpol_dt[13,:] = dpol_dt[13,:] .+ (aR_e .- aR_mu)
    dpol_dt[18,:] = dpol_dt[18,:] .+ (dpol_dt[2,:] *root3by2)

    dpol_dt[3:4,:] += 0.5*G2T5 .*(0.85*pol[3:4,:] .+ 0.75*pol[11:12,:])
    dpol_dt[6:7,:] += 0.5*G2T5 .*(3.56*pol[6:7,:])
    dpol_dt[8:9,:] += 0.5*G2T5 .*(2.5*pol[8:9,:])

    dpol_dt[11:12,:] += 0.5*G2T5 .*(0.85*pol[11:12,:] .+ 0.75*pol[3:4,:])
    dpol_dt[14:15,:] += 0.5*G2T5 .*(3.56*pol[14:15,:])
    dpol_dt[16:17,:] += 0.5*G2T5 .*(2.5*pol[16:17,:])

    # println(dpol_dt[10,:],dpol_dt[18,:])

    dpol_dt[:,:] = dpol_dt/(-T*Hubble(T)*(10^6))
end

###############################################################################   
Asym_e = 10^-10
Asym_mu = 10^-10

# pol_ini =  repeat([2*(2-(Asym_mu+Asym_e))/3,4/3,0,0,0,0,0,0,0,((2-(Asym_mu+Asym_e))/sqrt(3)),0,0,Asym_mu-Asym_e,0,0,0,0,(2/sqrt(3)),Asym_e/Lscale,Asym_mu/Lscale],1,lenx)
pol_ini =  repeat([4/3,2*(2-(Asym_mu+Asym_e))/3,0,0,0,0,0,0,0,(2/sqrt(3)),0,0,Asym_mu-Asym_e,0,0,0,0,((2-(Asym_mu+Asym_e))/sqrt(3)),Asym_e/Lscale,Asym_mu/Lscale],1,lenx)
# pol_ini =  repeat([4/3,4/3,0,0,0,0,0,0,0,(2/sqrt(3)),0,0,0,0,0,0,0,(2/sqrt(3)),Asym_e/Lscale,Asym_mu/Lscale],1,lenx)
# pol_ini =  repeat([0,0,0,0,0,0,0,0,0,0,0,0,Asym_mu-Asym_e,0,0,0,0,0,Asym_e/Lscale,Asym_mu/Lscale],1,lenx)
# pol_ini = convert(Matrix{Float16},pol_ini)
tspan = (60e6,0.1e6)

prob = ODEProblem(equation!,pol_ini,tspan,v)
# using Sundials
# sol = solve(prob,CVODE_BDF(),adaptive=false,dtmax=1e4)
# using ODEInterfaceDiffEq
# sol = solve(prob,dopri5(),adaptive=false,maxiters=1e6,dtmax=1e4)
@time sol = solve(prob, TRBDF2(),dtmax=1e3, maxiters=1e7 ,abstol=1e-7)#,reltol=1e-5)
# sol = solve(prob, IRKGL16(),saveat=np.linspace(20e6,1e6,1000),dtmax=1e3,dt=100.0,maxiters=1e6,abstol=1e-12,reltol=1e-12)

#TRBDF2

# using Sundials
# prob = DAEProblem(equation,pol_ini,tspan,p)
# sol = solve(prob,IDA())
# pol = zeros(Float64,20,lenx,size(sol.t)[1])
# for (i,_) in enumerate(sol.t)
#     pol[:,:,i] .= sol.u[i]
# end
np.save(name*"_T",sol.t)
np.save(name,sol);
println("done")

#plot(sol)
for i=1:lenx
    # plt.semilogx(sol.t,0.5*(sol[1,i,:] .+ sol[5,i,:] .+ (sol[10,i,:]/sqrt(3))))
    # plt.semilogx(sol.t,0.5*(sol[1,i,:] .- sol[5,i,:] .+ (sol[10,i,:]/sqrt(3))))
    # plt.semilogx(sol.t,0.5*(sol[2,i,:] .+ sol[13,i,:] .+ (sol[18,i,:]/sqrt(3))))
    # plt.semilogx(sol.t,0.5*(sol[2,i,:] .- sol[13,i,:] .+ (sol[18,i,:]/sqrt(3))))
    plt.semilogx(sol.t,0.5*(sol[1,i,:] .- (2*sol[10,i,:]/sqrt(3))))
    plt.semilogx(sol.t,0.5*(sol[2,i,:] .- (2*sol[18,i,:]/sqrt(3))))
end
plt.xlim(tspan[1],tspan[end])
display(gcf())
plt.close()
# for i=1:lenx
#     # plt.semilogy(sol.t,abs.(0.5*(pol[1,i,:] .+ pol[5,i,:] .+ (pol[10,i,:]/sqrt(3))) - 0.5*(pol[2,i,:] .+ pol[13,i,:] .+ (pol[18,i,:]/sqrt(3)))))
#     # plt.semilogy(sol.t,abs.(0.5*(pol[1,i,:] .- pol[5,i,:] .+ (pol[10,i,:]/sqrt(3))) - 0.5*(pol[2,i,:] .- pol[13,i,:] .+ (pol[18,i,:]/sqrt(3)))))
#     plt.semilogy(sol.t,abs.(sol[19,i,:])*1e-10)
#     plt.semilogy(sol.t,abs.(sol[20,i,:])*1e-10)
# end
# plt.xlim(tspan[1],tspan[end])
# display(gcf())
# plt.close()
# println("size:",size(sol.t))
#plt.plot(sol.t,0.5*(sol[1,:]+ sol[5,:] + (sol[10,:]/sqrt(3))))

#pol = np.swapaxes(pol,0,1)
u_e= 0.5*( sol[1,1,:] .+ sol[5,1,:] .+ (sol[10,1,:]/(np.sqrt(3))))
u_mu= 0.5*( sol[1,1,:] .- sol[5,1,:] .+ (sol[10,1,:]/(np.sqrt(3))))
u_s= 0.5*( sol[1,1,:] .- 2*(sol[10,1,:]/(np.sqrt(3))))
au_e= 0.5*( sol[2,1,:] .+ sol[13,1,:] .+ (sol[18,1,:]/(np.sqrt(3))))
au_mu= 0.5*( sol[2,1,:] .- sol[13,1,:] .+ (sol[18,1,:]/(np.sqrt(3))))
au_s= 0.5*( sol[2,1,:] .- 2*(sol[18,1,:]/(np.sqrt(3))))
n_eff_m = (1.5 * sol[1,1,:]) .- 2
n_eff_a = (1.5 * sol[2,1,:]) .- 2

fig1,ax1 =plt.subplots(figsize=(10,7))
ax1.loglog(sol.t/10^6,abs.(sol[19,1,:]*1e-10), label="L_e")
ax1.loglog(sol.t/10^6,abs.(sol[20,1,:]*1e-10), label="L_mu")
#ax1.set_title('3 Flavors with delta_CP = 0',fontsize=20,bbox=dict(facecolor='white', alpha=1))
ax1.set_ylabel("Lepton Asymmetries",fontsize=20)
ax1.set_xlabel("Temperature MeV", fontsize=20)
ax1.tick_params(axis="x", labelsize=0)
ax1.tick_params(axis="y", labelsize=20)
# ax1.set_yscale('symlog',linthresh=1e-29)
# ax1.set_ylim(1e-29,1e-1)
ax1.set_xlim(sol.t[1]/10^6,sol.t[end]/10^6)
# ax1.set_xlim(30,29)

ax1.legend(shadow=true,fontsize=20)
#fig1.savefig("2.png")
ax12 = fig1.add_axes((0.125,-0.175,0.775,0.3),sharex=ax1)
ax12.semilogx(sol.t/10^6,(n_eff_m .+ n_eff_a)/2,label="N_eff")
ax12.tick_params(axis="x",labelsize=15)
ax12.tick_params(axis="y",labelsize=15)
ax12.set_xlabel("Temperature (MeV)", fontsize=20)
ax12.set_ylabel("N_eff", fontsize=20)
# ax12.plot(sol.t/1e6,n_eff_a)
# ax12.set_ylim(bottom=-0.1,top= np.max((n_eff_m .+ n_eff_a)/2))
display(gcf())
plt.close()