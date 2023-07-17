using Plots
using LaTeXStrings
using LogarithmicNumbers #Not really needed

function cross_section(m, alpha)
    coupling_ratio = running_coupling_from_scale(2*m, Mtop, alpha_W_Mtop, 19/6-2/3*Ndark)/alpha #beta0 = 19/6 for the weak interaction below m_quark
    sigma0 = (NdarkAdjoint*(NdarkAdjoint-1)/(24*Ndark*Ndark) + coupling_ratio*(2*NdarkAdjoint/(3*Ndark) + 121/24*coupling_ratio))*pert_acs(alpha, m)/Ndark
end

function glueball_decay(m_glueball, m_quark, Alpha_DM)
    Alpha_DM_GB_decay = running_coupling_from_scale(m_glueball, m_quark, Alpha_DM, 11*Ndark/3) #Dark gauge coupling at the mass scale of the GBs
    Alpha_weak_GB_decay = running_coupling_from_scale(m_glueball, Mtop, alpha_W_Mtop, 19/6) #Weak gauge coupling at the mass scale of the GBs
    decay_const_GB = 3.06*m_glueball^3/(4*pi*Alpha_DM_GB_decay) #decay constant of the gluon after Juknevich. 
    Gamma_GB = (Alpha_weak_GB_decay*Alpha_DM_GB_decay)^2/(pi*m_quark^8)*1/1200*m_glueball^3*(decay_const_GB)^2 #Glueball decay rate after Juknevich. 
end

include("../FreezeOut.jl") #Include important functions regarding freeze-out
include("../SqueezeOut.jl") #Include important functions regarding seeze-out and GB decay

#Define constant physics parameters
const Ndark = 3. #Dark SU(N) parameter N
const NdarkAdjoint = Ndark*Ndark-1
const kappa = (Ndark*Ndark-1)*(Ndark*Ndark-2)/(16*Ndark^3)

#Running of SU(2)L gauge coupling
const MZ = 91.1876 # Z boson mass in GeV
const Mtop = 173 # Top quark mass in GeV
const alpha_W_MZ = 0.0342556 # Source?!
const alpha_W_Mtop = 0.0334 # The weak gauge coupling at the top quark mass (All from arxiv: 1307.3536)

#Constant parameters for entropy dilution
const R_max = 2.5E-4 #highest possible entropy ratio after the PT
const BBN_lifetime = 6.58*1E-25 #Lower bound on glueball decay rate.
const Relic_abundance_limit = 0.12 

const g_quark = 4*Ndark*3 #degeneracy of the Dirac quark: (Spin x Particle-Antiparticle) x DarkColour x weak multiplicity
gDM = 4*3 #For testing purposes

#Define physics parameters of the model (parameters in the loop)
Lambda_dQCD = 100
m_quark = Lambda_dQCD*40
Alpha_DM = running_coupling_from_pole(2*m_quark, Lambda_dQCD, 11*Ndark/3-2*3/3) #At the annihilation scale, the quark is active
Alpha_dark = running_coupling_from_pole(m_quark, Lambda_dQCD, 11*Ndark/3-2*3/3) #At the annihilation scale, the quark is active #dark gauge coupling at the mass scale m_quark
alpha = running_coupling_from_pole(2*m_quark, Lambda_dQCD, (11*Ndark-2)/3) #For testing purposes

BigConstant = bc_constant(m_quark)
sigma0 = cross_section(m_quark, Alpha_DM)
sigma_simple = kappa*pert_acs(alpha, m_quark) #For testing purposes
#Lambda_dQCD = Landau_pole(m_quark, Alpha_DM, 11) #beta0 = 11*Nc/3
Tcrit = 0.63*Lambda_dQCD #Temperature of the phase transition
x_PT = m_quark/Tcrit 

#Constants of FOPT (w/o factors of 1/Lambda for numerical convenience, cancel in the squeezeout step!)
R_pocket = pocket_radius(Lambda_dQCD)

#Define parameters of implicit Euler backward solution method
Delta_t = 1E-4
x_initial = 5
x_final = x_PT
#x_final = 1000
t_initial = log(x_initial)
t_final = log(x_final)

#First define initial conditions:
tvec = collect(t_initial:Delta_t:t_final)
xvec = exp.(tvec)
Npoints = length(tvec)

g_star_eff_vec = eff_dof_sqrt.(m_quark./xvec) #Effective degrees of freedom
g_effective = 106.75 #For testing purposes

EquilibriumYield = zeros(Npoints)
for i = 1:Npoints
    #EquilibriumYield[i] = Yeq(g_quark, h_eff_dof(m_quark/xvec[i]), xvec[i])
    EquilibriumYield[i] = Yeq(gDM, g_effective, xvec[i]) #For testing purposes
end

#Here the thermally averaged cross section is read in. 
#sigma_v_file = DataFrame(CSV.File("Sigma_eff.txt"))
#sigma_v_x_values = sigma_v_file[!,1]
#sigma_v_y_values = sigma_v_file[!,2]  #Not yet ready
#sigma_v_y_values = ones(500) #for evaluation w.o. SE and BSF

sigma_v_averaged = ones(Npoints) #Initialise array for interpolated annihaltion xs
#sigma_v_averaged_coeffs = sigma_v_interpolation(sigma0, sigma_v_x_values, sigma_v_y_values) #Cubic spline fit coefficient vector (beta, gamma, delta)

for i in 1:Npoints
    #sigma_v_averaged[i] = sigma_v_cspline(xvec[i], sigma_v_x_values, sigma_v_y_values, sigma0, sigma_v_averaged_coeffs)
    #sigma_v_averaged[i] = sigma0
    sigma_v_averaged[i] = sigma_simple #For testing purposes
end

Wx = zeros(Npoints)
Yx = zeros(Npoints)
Yx[1] = EquilibriumYield[1]
Wx[1] = log(Yx[1])

#Solution to the Boltzmann equation for the first freeze-out
for i = 2:Npoints
    W_old = Wx[i-1]
    #Wx[i] = Newton_Raphson_step(tvec[i], W_old, 0.5*BigConstant*g_star_eff_vec[i]*Delta_t*exp(-tvec[i])*sigma_v_averaged[i], g_quark, h_eff_dof(m_quark/xvec[i]))
    Wx[i] = Newton_Raphson_step(xvec[i], W_old, 0.5*BigConstant*g_effective*Delta_t/xvec[i]*sigma_v_averaged[i], gDM, g_effective) #For testing purposes
end
Yx = exp.(Wx)

#=
results_file = open("FO_data.csv", "w")
IOStream("FO_data.csv")
write(results_file, "x, Y \n")

for i in 1:1000:Npoints
    write(results_file, join((xvec[i],Yx[i]), ","),"\n")
end
close(results_file)
=#
### FOPT: Squeezeout step ###
Yx_squeezeout = 1.5/pi*sqrt(15*Yx[Npoints]/(2*pi*h_eff_dof(Tcrit)*R_pocket^3))

if m_quark/Lambda_dQCD < 30 #Intermediate and strong coupling regime
    ### Include a new 'freeze-out' solution stage, now for baryons
    x_initial2 = x_PT
    x_final2 = 1000
    t_initial2 = log(x_initial2)
    t_final2 = log(x_final2)

    #First define initial conditions:
    tvec2 = collect(t_initial2:Delta_t:t_final2)
    xvec2 = exp.(tvec2)
    Npoints2 = length(tvec2)

    EquilibriumYield2 = zeros(Npoints2)
    for i = 1:Npoints2
        #EquilibriumYield[i] = Yeq(g_quark, h_eff_dof(m_quark/xvec[i]), xvec[i])
        EquilibriumYield2[i] = Yeq(gDM, g_effective, xvec2[i]) #For testing purposes
    end

    sigma_v_averaged2 = ones(Npoints2) #Initialise array for interpolated annihaltion xs
    #sigma_v_averaged_coeffs = sigma_v_interpolation(sigma0, sigma_v_x_values, sigma_v_y_values) #Cubic spline fit coefficient vector (beta, gamma, delta)

    for i in 1:Npoints2
        #sigma_v_averaged[i] = sigma_v_cspline(xvec[i], sigma_v_x_values, sigma_v_y_values, sigma0, sigma_v_averaged_coeffs)
        #sigma_v_averaged[i] = sigma0
        sigma_v_averaged2[i] = sigma_simple #For testing purposes
    end

    Wx2 = zeros(Npoints2)
    Yx2 = zeros(Npoints2)
    Yx2[1] = Yx_squeezeout
    Wx2[1] = log(Yx2[1])

    #Solution to the Boltzmann equation for the second "freeze-out"
    for i = 2:Npoints2
        W_old = Wx2[i-1]
        #Wx[i] = Newton_Raphson_step(tvec[i], W_old, 0.5*BigConstant*g_star_eff_vec[i]*Delta_t*exp(-tvec[i])*sigma_v_averaged[i], g_quark, h_eff_dof(m_quark/xvec[i]))
        Wx2[i] = Newton_Raphson_step(xvec2[i], W_old, 0.5*BigConstant*g_effective*Delta_t/xvec2[i]*sigma_v_averaged2[i], gDM, g_effective) #For testing purposes
    end
    Yx2 = exp.(Wx2)
end

### Entropy dilution due to glueball decay ###
m_glueball = 7*Lambda_dQCD #Mass of the lightest 0++ glueball
x_freeze_out = GB_freeze_out_estimate(1, m_glueball, R_max) #freeze-out of dark gluons. xft = 1 is the initial guess 

#Matter-radiation equality
Y_GB = R_max/x_freeze_out #Relic yield of dark gluons
T_MR = 4/3*m_glueball*Y_GB # Matter-radiation equality temperature, after which GB dominate the energy content
x_MR = m_quark/T_MR

#Glueball decay
Gamma_GB = glueball_decay(m_glueball, m_quark, Alpha_DM)

#Entropy dilution
dil_fac = (1 + 1.65*g_average(1e-5, T_MR, 10)*cbrt(T_MR^4/(Gamma_GB*Mpl))^2)^(-0.75)
if m_quark/Lambda_dQCD < 30 #Intermediate and strong coupling regime
    Yx_dilution = dil_fac*Yx2[Npoints2]
else #Coulomb regime
    Yx_dilution = dil_fac*Yx_squeezeout
end
relic_abundance_Yield_limit = Relic_abundance_limit*rho_crit/(reduced_Hubble_squared*s0*m_quark)

plot(xvec, EquilibriumYield, title="V model: m = $(m_quark/1000) TeV, " * L"\Lambda = " * "$Lambda_dQCD GeV",
    label=[L"Y_{eq}(x)" L"Y_{f.o.}(x)"], minorticks = 10, minorgrid = true, xlabel="x = m/T", ylabel="Y(x)",
    xaxis=:log, yaxis=:log, xlims = (x_initial,1E4), ylims = (1E-23, 1E-1), linewidth = 3)
plot!(xvec, Yx,  label=L"Y_{quark}(x)", linewidth = 3)

if m_quark/Lambda_dQCD < 30 #Intermediate and strong coupling regime
    plot!(xvec2, Yx2,  label=L"Y_{baryon}(x)", linewidth = 3)
else #Coulomb regime
    plot!([x_PT; 1E3],[Yx_squeezeout; Yx_squeezeout], label = L"Y_{s.q.}", linewidth = 3)
end
plot!([1E3; 1E4],[Yx_dilution; Yx_dilution], label = L"Y_{dil}", linewidth = 3)
hline!([relic_abundance_Yield_limit], linestyle=:dash, label = L"\Omega h^2 = 0.12", linewidth = 3)
savefig("ThermalHistory_conf_AFTER_fo.png")