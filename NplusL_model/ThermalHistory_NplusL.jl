include("NplusL_Header.jl") #Include important functions 

Lambda_dQCD = 100
m_N = 100*Lambda_dQCD
mass_delta = 4
m_L = (mass_delta+1)*m_N 
ydark = 0.1 #dark Yukawa coupling. No running implemented. 

Alpha_DM = running_coupling_from_pole(m_N, Lambda_dQCD, 11*Ndark/3)

Tcrit = 1.91*Lambda_dQCD #Temperature of the phase transition
x_PT = m_N/Tcrit 

#Constants of FOPT (w/o factors of 1/Lambda for numerical convenience, cancel in the squeezeout step!)
R_pocket = pocket_radius(Lambda_dQCD)


###Freeze-Out
#Define parameters of implicit Euler backward solution method
Delta_t = 1E-4
x_initial = 1E-1
x_final = x_PT
t_initial = log(x_initial)
t_final = log(x_final)

#First define initial conditions:
tvec = collect(t_initial:Delta_t:t_final)
xvec = exp.(tvec)
Npoints = length(tvec)

g_star_eff_vec = eff_dof_sqrt.(m_N./xvec) #Effective degrees of freedom
h_eff_dof_vec = h_eff_dof.(m_N./xvec) #Effective entropic dofs
BigConstant = bc_constant(m_N)

#Set Equilibrium yield (for initial condition and for plotting)
EquilibriumYieldDM = zeros(Npoints)
for i in 1:Npoints
    EquilibriumYieldDM[i] = Yeq_DM_coannihilation(g_N, g_L, h_eff_dof_vec[i], m_L/m_N, xvec[i])
end

#Define the initial conditions
Wx = zeros(Npoints)
Yx = zeros(Npoints)
Yx[1] = EquilibriumYieldDM[1]
Wx[1] = log(Yx[1])

mod_BC = 0.5*BigConstant*Delta_t #modified Big Constant
sigma_v_averaged = zeros(Npoints) #Calculation of the averaged cross section
(sigma_NN, sigma_NL, sigma_LL) = cross_section_contributions(m_N, m_L, Lambda_dQCD, ydark)

#BEGIN BLOCK Interpolation
num_raw_data = 100
sigma_v_raw_data = ones(num_raw_data)
x_raw_data_vec = zeros(num_raw_data)
x_raw_data_vec[1] = xvec[1]
sigma_v_raw_data_first_entry = cross_section(xvec[1], m_N, m_L, Lambda_dQCD, ydark, sigma_NN, sigma_NL, sigma_LL)
for i in 2:num_raw_data
    raw_ind = floor(Int, Npoints/num_raw_data*(i-1) + 1)
    x_raw_data_vec[i] = xvec[raw_ind]
    sigma_v_raw_data[i] = cross_section(xvec[raw_ind], m_N, m_L, Lambda_dQCD, ydark, sigma_NN, sigma_NL, sigma_LL)/sigma_v_raw_data_first_entry
end

interpolation_coeffs = sigma_v_interpolation(x_raw_data_vec, sigma_v_raw_data) ###array of interpolation coefficients
####

for i in 1:Npoints
    sigma_v_averaged[i] = sigma_v_cspline(xvec[i], x_raw_data_vec, sigma_v_raw_data, sigma_v_raw_data_first_entry, interpolation_coeffs)
    #sigma_v_averaged[i] = cross_section(xvec[i], m_N, m_L, Lambda_dQCD, ydark, sigma_NN, sigma_NL, sigma_LL)
end

#Solution to the Boltzmann equation for the first freeze-out
for i = 2:Npoints
    W_old = Wx[i-1]
    Wx[i] = Newton_Raphson_step_coannihilation(xvec[i], W_old, mod_BC*sigma_v_averaged[i]/xvec[i], g_N, g_L, h_eff_dof_vec[i], m_L/m_N)
end
Yx = exp.(Wx)
Yfo = Yx[Npoints]

#Yfo = coannihilation_freeze_out(x_PT, m_N, m_L, Lambda_dQCD, ydark)

### FOPT: Squeezeout step ###
Yx_squeezeout = 1.5/pi*sqrt(15*Yfo/(2*pi*h_eff_dof(Tcrit)*R_pocket^3))

### Entropy dilution due to glueball decay ###
m_glueball = 7*Lambda_dQCD #Mass of the lightest 0++ glueball
x_freeze_out = GB_freeze_out_estimate(1, m_glueball, R_max) #freeze-out of dark gluons

#Matter-radiation equality
Y_GB = R_max/x_freeze_out #Relic yield of dark gluons
T_MR = 4/3*m_glueball*Y_GB # Matter-radiation equality temperature, after which GB dominate the energy content

#GB decay
Gamma_GB = glueball_decay(m_glueball, Lambda_dQCD, m_N, m_L, ydark)

#Entropy dilution
dil_fac = (1 + 1.65*g_average(1e-5, T_MR, 10)*cbrt(T_MR^4/(Gamma_GB*Mpl))^2)^(-0.75)
#x_dilution = x_PT*100
Yx_dilution = dil_fac*Yx_squeezeout

### Calculation of the relic abundance ###
x_today = m_N/T0
Omega_relic = reduced_Hubble_squared*Yx_dilution*s0*m_N/rho_crit
relic_abundance_Yield_limit = Relic_abundance_limit*rho_crit/(reduced_Hubble_squared*s0*m_N)

plot(xvec, sigma_v_averaged, title="Effective cross section", label=L"\sigma v"*" interpolated", minorticks = 10,
 minorgrid = true, xlabel = L"x = m_N/T", ylabel=L"\sigma v", xaxis=:log, yaxis=:log, linewidth = 3)
scatter!(x_raw_data_vec, sigma_v_raw_data_first_entry*sigma_v_raw_data, label=L"\sigma v"*" measured")
savefig("Interpolated_cross_section.png")

plot(xvec, EquilibriumYieldDM, title="N+L: m = $(m_N/1000) TeV, " * L"\Lambda = " * "$Lambda_dQCD GeV, "
 * L"y_d = " * "$ydark, " * L"\delta = " * "$mass_delta",
    label=[L"Y_{eq}(x)"], minorticks = 10, minorgrid = true, xlabel=L"x = m_N/T", ylabel="Y(x)",
    xaxis=:log, yaxis=:log, xlims = (x_initial,1E4), ylims = (1E-23, 1E-1), linewidth = 3)
plot!(xvec, Yx,  label=L"Y_{quark}(x)", linewidth = 3)
plot!([x_PT; 1E3],[Yx_squeezeout; Yx_squeezeout], label = L"Y_{s.q.}", linewidth = 3)
#plot!(xvec2, Yx2,  label=L"Y_{baryon}(x)", linewidth = 3)
plot!([1E3; 1E4],[Yx_dilution; Yx_dilution], label = L"Y_{dil}", linewidth = 3)
hline!([relic_abundance_Yield_limit], linestyle=:dash, label = L"\Omega h^2 = 0.12", linewidth = 3)
savefig("ThermalHistory_NplusL_conf_AFTER_fo.png")
