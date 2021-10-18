using PyCall
sp = pyimport("sympy")
a = Array([[1 0 0];[0 2 0];[0 0 3]])
b = Array([[2 2 4];[3 2 2];[3 6 3]])
c = ones(3,3)
(a^2) * b

function collide(pol)
    G2T5 = (4.26984*10^-17) * T^5
    R_e = 0.5*G2T5*(13.24*(Nu_e -1) + (aNu_e -1))
    R_mu = 0.5*G2T5*(9.44*(Nu_mu -1) + 0.56*(aNu_mu -1))
    aR_e = 0.5*G2T5*(13.24*(aNu_e -1) + (Nu_e -1))
    aR_mu = 0.5*G2T5*(9.44*(aNu_mu -1) + 0.56*(Nu_mu -1))

    dpol_dt[1,:] = (2/3)*(R_e .+ R_mu)
    dpol_dt[5,:] = dpol_dt[5,:] .+ (R_e .- R_mu)
    dpol_dt[10,:] = dpol_dt[10,:] .+ (dpol_dt[1,:] *root3by2)

    dpol_dt[2,:] = (2/3)*(aR_e .+ aR_mu)
    dpol_dt[13,:] = dpol_dt[13,:] .+ (aR_e .- aR_mu)
    dpol_dt[18,:] = dpol_dt[18,:] .+ (dpol_dt[2,:] *root3by2)

    dpol[3:4,:] += 0.5*G2T5*(11.28*pol[3:4,:] + 0.75*pol[11:12,:])
    dpol[6:7,:] += 0.5*G2T5*(3.56*pol[6:7,:])
    dpol[8:9,:] += 0.5*G2T5*(2.5*pol[8:9,:])

    dpol[11:12,:] += 0.5*G2T5*(11.28*pol[11:12,:] + 0.75*pol[3:4,:])
    dpol[14:15,:] += 0.5*G2T5*(3.56*pol[14:15,:])
    dpol[16:17,:] += 0.5*G2T5*(2.5*pol[16:17,:])