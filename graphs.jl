pol_t =pol ; T = sol.t
Nu_e= 0.5*(pol_t[1,:,:] .+ pol_t[5,:,:] .+ (pol_t[10,:,:]/sqrt3))
aNu_e= 0.5*(pol_t[2,:,:] .+ pol_t[13,:,:] .+ (pol_t[18,:,:]/sqrt3))
Nu_mu= 0.5*(pol_t[1,:,:] .- pol_t[5,:,:] .+ (pol_t[10,:,:]/sqrt3))
aNu_mu= 0.5*(pol_t[2,:,:] .- pol_t[13,:,:] .+ (pol_t[18,:,:]/sqrt3))
Nu_s= 0.5*(pol_t[1,:,:] .- 2*(pol_t[10,:,:]/sqrt3))
aNu_s= 0.5*(pol_t[2,:,:] .- 2*(pol_t[18,:,:]/sqrt3))

# Nu_mu_norm, aNu_mu_norm = take_norm.([Nu_mu,aNu_mu])
f_veL = zeros(lenx,size(T)[1])
f_vmuL = zeros(lenx,size(T)[1])
f_ve_L = zeros(lenx,size(T)[1])
f_vmu_L = zeros(lenx,size(T)[1])
dpol_dt =zeros(20,lenx,size(T)[1])
v= zeros(16,lenx,size(T)[1])
pot_allv= zeros(16,lenx)
# Nu_e=zeros(lenx,size(sol.t)[1])
# Nu_mu=zeros(lenx,size(sol.t)[1])
# aNu_e=zeros(lenx,size(sol.t)[1])
# aNu_mu=zeros(lenx,size(sol.t)[1])
for (index,j) in enumerate(T)
    # Nu_e[:,index],Nu_mu[:,index],aNu_e[:,index],aNu_mu[:,index] = equation(dpol_dt[:,:,index],sol[:,:,index],pot_allv,j)
    # dpol_dt[:,:,index] = equation(dpol_dt[:,:,index],pol_t[:,:,index],pot_allv,j)
    v[:,:,index] = repeat(vacuum_potential(x,j/(10^6)),2)
    # Nu_mu_norm, aNu_mu_norm = take_norm.([Nu_mu[:,index],aNu_mu[:,index]])
    f_veL[:,index],f_ve_L[:,index] = forward_v4e(j/(10^6),Nu_e[:,index],aNu_e[:,index],0)
    f_vmuL[:,index],f_vmu_L[:,index] =forward_v(j/(10^6),Nu_mu[:,index],aNu_mu[:,index],0)
    v[3,:,index] .+= (f_veL[:,index] .- f_vmuL[:,index])/2
    v[8,:,index] .+= (f_veL[:,index] .+ f_vmuL[:,index])/twoxsqrt3
end

# for (index,j) in enumerate(sol.t)
#     # plt.plot(x,Nu_mu[:,index])
#     # plt.plot(x,Nu_e[:,index])
#     # plt.plot(x,Nu_s[:,index])
#     # v[3,:,index]= v[3,:,index] .+ (- f_vmuL[:,index])/2
#     # v[8,:,index]= v[8,:,index] .+ (f_vmuL[:,index])/twoxsqrt3
# end
# plt.xlim(x[1],50)
# display(gcf())
# plt.close()

for i=1:lenx
    # plt.plot(T,Nu_e[i,:])
    # plt.plot(T,aNu_e[i,:])
    # plt.plot(T,Nu_mu[i,:])
    # plt.plot(T,aNu_mu[i,:])
    # plt.semilogy(T,abs.(f_veL[i,:]))
    # plt.semilogy(T,abs.(f_vmuL[i,:]))
    # plt.semilogy(sol.t,abs.(sol[8,i,:]))
    # plt.semilogy(sol.t,abs.(sol[9,i,:]))
    # plt.semilogy(sol.t,abs.(sol[10,i,:]))
    # plt.semilogy(sol.t,abs.(dpol_dt[5,i,:]),label="5")
    # plt.semilogy(sol.t,abs.(dpol_dt[6,i,:]),label="6")
    # plt.semilogy(sol.t,abs.(dpol_dt[7,i,:]),label="7")
    # plt.semilogy(sol.t,abs.(dpol_dt[8,i,:]),label="8")
    # plt.semilogy(sol.t,abs.(dpol_dt[9,i,:]),label="9")
    # plt.semilogy(sol.t,abs.(dpol_dt[10,i,:]),label="10")
    plt.semilogy(T,abs.(v[3,i,:] .+ v[8,i,:]*sqrt3))
    plt.semilogy(T,abs.(-v[3,i,:] .+ v[8,i,:]*sqrt3))
    # plt.plot(sol.t,0.5 .*(sol[2,i,:]- sol[13,i,:] + (sol[18,i,:] ./sqrt(3))))
    # plt.plot(sol.t,0.5 .*(sol[1,i,:] - 2 .*(sol[10,i,:] ./sqrt(3))))
    # plt.plot(sol.t,0.5 .*(sol[2,i,:] - 2 .*(sol[18,i,:] ./sqrt(3))))
end
plt.legend()
plt.xlim(tspan[1],tspan[end])
# plt.xlim(7e6,4.5e6)
# plt.ylim(0,2)
display(gcf())
plt.close()
# typeof(sol.t)
# pol_t= permutedims(np.load("august_track.npy"),[2,1,3])
# dpol_dt[:,:,50] = equation(dpol_dt[:,:,50],pol_t[:,:,50],pot_allv,sol.t[50])


# function equation(pol::Matrix{Z},v,T) where {Z <: Number}
#     T=T/10^6
#     # pol = Vector{Z}(20,lenx) 
#     # cross1 = Vector{Z}(1,lenx)
#     #v=zeros(Real,18,lenx)
#     # pol = reshape(pol_in,20,lenx)
#     dpol_dt = zeros(20,lenx)

#     Nu_e= 0.5*(pol[1,:] .+ pol[5,:] .+ (pol[10,:]/sqrt3))
#     Nu_mu= 0.5*(pol[1,:] .- pol[5,:] .+ (pol[10,:]/sqrt3))

#     aNu_e= 0.5*(pol[2,:] .+ pol[13,:] .+ (pol[18,:]/sqrt3))
#     aNu_mu= 0.5*(pol[2,:] .- pol[13,:] .+ (pol[18,:]/sqrt3))

#     Nu_e_norm, aNu_e_norm, Nu_mu_norm, aNu_mu_norm = take_norm.([Nu_e,aNu_e,Nu_mu,aNu_mu])

#     L_e = (pol[19,1]).*Lscale .*0
#     L_mu = (pol[20,1]).*Lscale .*0 
#     chem_p_e = -3.628*pi*sinh.(asin.(-1.2087*(pol[19,1] .*Lscale))/3) .*0
#     chem_p_mu= -3.628*pi*sinh.(asin.(-1.2087*(pol[20,1] .*Lscale))/3) .*0

#     v[:,:] = repeat(vacuum_potential(x,T),2)

#     G2T5 = GGG*(T^5)

#     f_veL,f_ve_L = forward_v4e(T,Nu_e_norm,aNu_e_norm,L_e)
#     f_vmuL,f_vmu_L = forward_v(T,Nu_mu_norm,aNu_mu_norm,L_mu)

#     cross_real, cross_imag = (1.49851*(10^-6)*(T^3)).*[pol[3,:] .- pol[11,:], pol[4,:] .- pol[12,:]] 

#     # println(typeof(v))
#     v[1,:]= v[1,:] .+ cross_real
#     v[2,:]= v[2,:] .- cross_imag
#     v[3,:]= v[3,:] .+ (f_veL .- f_vmuL)/2
#     v[8,:]= v[8,:] .+ (f_veL .+ f_vmuL)/twoxsqrt3

#     v[9,:] = v[9,:] .- cross_real
#     v[10,:]= v[10,:] .+ cross_imag
#     v[11,:]= v[11,:] .+ (f_ve_L .- f_vmu_L)/2
#     v[16,:]= v[16,:] .+ (f_ve_L .+ f_vmu_L)/twoxsqrt3

#     # println(v)
#     # open("file.txt","a") do io
#     #     println(io,"a=",v[3,:])
#     #     println(io,"b=",v[8,:])
#     # end
#     dpol_dt[3:18,:] = GMcross(2*v,pol)

#     R_e = G2T5.*(C_e.*(fermi(x,chem_p_e)./f0 .- Nu_e))
#     R_mu = G2T5.*(C_mu.*(fermi(x,chem_p_mu)./f0 .- Nu_mu))
#     aR_e = G2T5.*(C_e.*(fermi(x,-chem_p_e)./f0 .- aNu_e))
#     aR_mu = G2T5.*(C_mu.*(fermi(x,-chem_p_mu)./f0 .- aNu_mu))  

#     dpol_dt[19,:] .= 0#trapz(x,((dpol_dt[5,:].-dpol_dt[13,:]) .+ ((dpol_dt[10,:].-dpol_dt[18,:]).*sqrt3)) .*f0x2)/(16*zeta3*Lscale)
#     dpol_dt[20,:] .= 0#trapz(x,((dpol_dt[13,:].-dpol_dt[5,:]) .+ ((dpol_dt[10,:].-dpol_dt[18,:]).*sqrt3)) .*f0x2)/(16*zeta3*Lscale)

#     dpol_dt[1,:] = (2/3)*(R_e .+ R_mu)
#     dpol_dt[5,:] += (R_e .- R_mu)
#     dpol_dt[10,:] += dpol_dt[1,:]*root3by2

#     dpol_dt[6:7,:] .= dpol_dt[6:7,:] .- (coeffs[1,:][:,[CartesianIndex()]] .*G2T5[[CartesianIndex()],:]) .* pol[6:7,:]
#     dpol_dt[8:9,:] .= dpol_dt[8:9,:] .- (coeffs[2,:][:,[CartesianIndex()]] .*G2T5[[CartesianIndex()],:]) .* pol[8:9,:]
#     dpol_dt[3:4,:] .= dpol_dt[3:4,:] .-(coeffs[3,:][:,[CartesianIndex()]] .*G2T5[[CartesianIndex()],:]) .* pol[3:4,:]
#     dpol_dt[14:15,:] .= dpol_dt[14:15,:] .- (coeffs[4,:][:,[CartesianIndex()]] .*G2T5[[CartesianIndex()],:]) .* pol[14:15,:]
#     dpol_dt[16:17,:] .= dpol_dt[16:17,:] .- (coeffs[5,:][:,[CartesianIndex()]] .*G2T5[[CartesianIndex()],:]) .* pol[16:17,:]
#     dpol_dt[11:12,:] .= dpol_dt[11:12,:] .- (coeffs[6,:][:,[CartesianIndex()]] .*G2T5[[CartesianIndex()],:]) .* pol[11:12,:]

#     dpol_dt[3:4,:] .= 0#dpol_dt[3:4,:] .- (C_dash[1,:][:,[CartesianIndex()]] .*G2T5[[CartesianIndex()],:]) .* pol[11:12,:]
#     dpol_dt[11:12,:] .= 0#dpol_dt[11:12,:] .- (C_dash[2,:][:,[CartesianIndex()]] .*G2T5[[CartesianIndex()],:]) .* pol[3:4,:]
#     #k = [dpol_dt[6:7,:],dpol_dt[8:9,:],dpol_dt[3:4,:],dpol_dt[14:15,:],dpol_dt[16:17,:],dpol_dt[11:12,:]]

#     #dpol_dt[6:7,:],dpol_dt[8:9,:],dpol_dt[3:4,:],dpol_dt[14:15,:],dpol_dt[16:17,:],dpol_dt[11:12,:] = k .- (coeffs[:,[CartesianIndex()]] .*G2T5[[CartesianIndex()],:])[:,[CartesianIndex()],:].* permutedims(cat(pol[6:7,:],pol[8:9,:],pol[3:4,:],pol[14:15,:],pol[16:17,:],pol[11:12,:],dims=3),(3,1,2))
#     #dpol_dt[3:4,:],dpol_dt[11:12,:] = (dpol_dt[3:4,:],dpol_dt[11:12,:]) .- (C_dash[:,[CartesianIndex()]] .* G2T5[[CartesianIndex()],:])[:,[CartesianIndex()],:].* permutedims(cat(pol[11:12,:],pol[3:4,:],dims=3),(3,1,2))

#     dpol_dt[2,:] = (2/3)*(aR_e .+ aR_mu)
#     dpol_dt[13,:] += (aR_e .- aR_mu)
#     dpol_dt[18,:] += dpol_dt[2,:]*root3by2
#     #dpol_dt[:,2:4],dpol_dt[:,10:12] = np.array([dpol_dt[:,2:4],dpol_dt[:,10:12]])/1e-16

#     dpol_dt = dpol_dt ./(-Hubble(T)*(T*10^6))
#     return dpol_dt
# end
