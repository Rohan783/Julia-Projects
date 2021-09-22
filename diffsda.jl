# using ForwardDiff
# f(x::Float64) = x^2
# function f_diff(x) 
#     ForwardDiff.derivative(f, x)
# end
# f_diff(2.0)

function test(a,f,b,sub)
    e = (a[:,:,[CartesianIndex()]] .*f[[CartesianIndex()],[CartesianIndex()],:]) .* permutedims(cat(b,b,b,b,b,b,dims=3),(3,1,2))
    sub[1] = e[1,:,:]
    sub[2] = e[2,:,:]
    sub[3] = e[3,:,:]
    sub[4] = e[4,:,:]
    sub[5] = e[5,:,:]
    sub[6] = e[6,:,:]
    return sub
    #return [e[i,:,:] for i=1:size(e)[1]]
end


function test2(a,f,b,sub)
    # e = (a[:,:,[CartesianIndex()]] .*f[[CartesianIndex()],[CartesianIndex()],:]) .* permutedims(cat(b,b,b,b,b,b,dims=3),(3,1,2))
    sub[1] = (a[1,:][:,[CartesianIndex()]]*f[[CartesianIndex()],:]) .* b
    sub[2] = (a[2,:][:,[CartesianIndex()]]*f[[CartesianIndex()],:]) .* b
    sub[3] = (a[3,:][:,[CartesianIndex()]]*f[[CartesianIndex()],:]) .* b
    sub[4] = (a[4,:][:,[CartesianIndex()]]*f[[CartesianIndex()],:]) .* b
    sub[5] = (a[5,:][:,[CartesianIndex()]]*f[[CartesianIndex()],:]) .* b
    sub[6] = (a[6,:][:,[CartesianIndex()]]*f[[CartesianIndex()],:]) .* b
    return sub
    #return [e[i,:,:] for i=1:size(e)[1]]
end

a=reshape(collect(1:12),6,2)
b= ones(2,5)
f= ones(5)*2
ddol =Array{Any}(undef,6)
#d = repeat(coeffs,1,2)[:,:,[CartesianIndex()]] .*GGG[[CartesianIndex()],[CartesianIndex()],:]
c= (a[:,:,[CartesianIndex()]] .*f[[CartesianIndex()],[CartesianIndex()],:]) .* permutedims(cat(b,b,b,b,b,b,dims=3),(3,1,2))
@time test(a,f,b,ddol)
@time test2(a,f,b,ddol)