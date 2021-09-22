# pol_t= permutedims(np.load("august_track.npy"),[2,1,3])
# dpol_dt[:,:,50] = equation(dpol_dt[:,:,50],pol_t[:,:,50],pot_allv,sol.t[50])

T=sol.t
# T=np.linspace(20e6,1e6,1000)
dpol_dt =zeros(20,lenx,size(T)[1])
v= zeros(16,lenx,size(T)[1])
pot_allv= zeros(16,lenx)
a = zeros(lenx,size(T)[1])
b = zeros(lenx,size(T)[1])
c = zeros(lenx,size(T)[1])
for (index,j) in enumerate(T[2:end])
    # Nu_e[:,index],Nu_mu[:,index],aNu_e[:,index],aNu_mu[:,index] = equation(dpol_dt[:,:,index],sol[:,:,index],pot_allv,j)
    dpol_dt[:,:,index] = equation!(dpol_dt[:,:,index],pol_t[:,:,index],pot_allv,j)
    # a[:,index+1],b[:,index+1],c[:,index+1] = equation(dpol_dt[:,:,index],pol_t[:,:,index],pot_allv,j)
end

for i=1:lenx
    # plt.plot(sol.t,Nu_e[i,:])
    # plt.plot(sol.t,aNu_e[i,:])
    # plt.plot(sol.t,Nu_mu[i,:])
    # plt.plot(sol.t,aNu_mu[i,:])
    # plt.semilogy(sol.t,abs.(f_vmuL[i,:]))
    # plt.semilogy(T,abs.(a[i,:]))
    # plt.semilogy(T,abs.(b[i,:]))
    # plt.semilogy(T,abs.(c[i,:]))
    # plt.semilogy(T,abs.(a[i,:]+b[i,:]))
    # plt.semilogy(sol.t,abs.(sol[10,i,:]))
    plt.semilogy(T,abs.(dpol_dt[5,i,:]),label="5")
    plt.semilogy(T,abs.(dpol_dt[6,i,:]),label="6")
    plt.semilogy(T,abs.(dpol_dt[7,i,:]),label="7")
    plt.semilogy(T,abs.(dpol_dt[8,i,:]),label="8")
    plt.semilogy(T,abs.(dpol_dt[9,i,:]),label="9")
    plt.semilogy(T,abs.(dpol_dt[10,i,:]),label="10")
    # plt.semilogy(sol.t,abs.(v[3,i,:] .+v[8,i,:]*sqrt3))
    # plt.semilogy(sol.t,abs.(-v[3,i,:] .+v[8,i,:]*sqrt3))
    # plt.plot(sol.t,0.5 .*(sol[2,i,:]- sol[13,i,:] + (sol[18,i,:] ./sqrt(3))))
    # plt.plot(sol.t,0.5 .*(sol[1,i,:] - 2 .*(sol[10,i,:] ./sqrt(3))))
    # plt.plot(sol.t,0.5 .*(sol[2,i,:] - 2 .*(sol[18,i,:] ./sqrt(3))))
end
plt.legend()
# plt.xlim(tspan[1],tspan[end])
# plt.yscale("symlog",linthresh=1e-13)
# plt.xlim(10e6,tspan[end])
plt.xlim(20e6,1e6)
# plt.ylim(0.8e-13,0.3e-9)
display(gcf())
plt.close()