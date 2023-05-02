using Plots

include("../FreezeOut.jl")

g_quark = 4 #degeneracy of the Dirac quark

m_quark = 10000
Alpha_DM = 0.1

BigConstant = bc_constant(m_quark)
sigma0 = pert_acs(Alpha_DM, m_quark)

#Define parameters of implicit Euler backward solution method
x_final = 1000

Delta_t = 1E-4
x_initial = 1E-1
t_initial = log(x_initial)
t_final = log(x_final)

#First define initial conditions:
tvec = collect(t_initial:Delta_t:t_final)
xvec = exp.(tvec)
Npoints = length(tvec)

g_star_eff_vec = eff_dof_sqrt.(m_quark./xvec) #Effective degrees of freedom

EquilibriumYield = zeros(Npoints)
for i = 1:Npoints
    EquilibriumYield[i] = Yeq(g_quark, h_eff_dof(m_quark/xvec[i]), xvec[i])
end

sigma_v_averaged = ones(Npoints) #Initialise array for interpolated annihaltion xs
#sigma_v_averaged_coeffs = sigma_v_interpolation(sigma0, sigma_v_x_values, sigma_v_y_values) #Cubic spline fit coefficient vector (beta, gamma, delta)

for i in 1:Npoints
    #sigma_v_averaged[i] = sigma_v_cspline(xvec[i], sigma_v_x_values, sigma_v_y_values, sigma0, sigma_v_averaged_coeffs)
    sigma_v_averaged[i] = sigma0
end

Wx = zeros(Npoints)
Yx = zeros(Npoints)
Yx[1] = EquilibriumYield[1]
Wx[1] = log(Yx[1])

#Solution to the Boltzmann equation for the first freeze-out
for i = 2:Npoints
    W_old = Wx[i-1]
    Wx[i] = Newton_Raphson_step(tvec[i], W_old, 0.5*BigConstant*Delta_t*sigma_v_averaged[i], g_quark, h_eff_dof(m_quark/xvec[i]))
end

Yx = exp.(Wx)

plot(xvec, Yx, xscale = :log10, yscale = :log10, xlabel = "x = m/T", ylabel = "Y(x)", title = "Freezeout", minorgrid = true, minorticks = 10)
savefig("Test.png")