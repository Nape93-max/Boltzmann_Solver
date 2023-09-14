using Plots
using LaTeXStrings

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
const g_Baryon = 4*6*Ndark #degeneracy of the baryon. 6 = 3*2*1 SU(2)L dofs and 4*Ndark is derived in the Overleaf document

g_effective = 106.75 #For testing purposes

#Define physics parameters of the model (parameters in the loop)
Lambda_dQCD_1 = 50
m_quark_1 = Lambda_dQCD_1*5
Lambda_dQCD_2 = 50
m_quark_2 = Lambda_dQCD_2*1E4
alpha_1 = running_coupling_from_pole(2*m_quark_1, Lambda_dQCD_1, (11*Ndark-2)/3) #For testing purposes
alpha_2 = running_coupling_from_pole(2*m_quark_2, Lambda_dQCD_2, (11*Ndark-2)/3) #For testing purposes

BigConstant_1 = bc_constant(m_quark_1)
BigConstant_2 = bc_constant(m_quark_2)
sigma_1 = pert_acs(alpha_1, m_quark_1)
sigma_12 = 4*pi/(alpha_1*m_quark_1*m_quark_1)
sigma_2 = cross_section(m_quark_2, alpha_2)
Tcrit_1 = 1.2*Lambda_dQCD_1 #Temperature of the phase transition
Tcrit_2 = 1.2*Lambda_dQCD_2 #Temperature of the phase transition
x_PT_1 = m_quark_1/Tcrit_1 
x_PT_2 = m_quark_2/Tcrit_2 

#Constants of FOPT (w/o factors of 1/Lambda for numerical convenience, cancel in the squeezeout step!)
R_pocket = pocket_radius(Lambda_dQCD_2)

#Define parameters of implicit Euler backward solution method
Delta_t = 1E-4
x_initial_1 = 0.5*m_quark_1*Lambda_dQCD_2/m_quark_2/Lambda_dQCD_1
x_initial_2 = 0.5
x_final_1 = x_PT_1
x_final_2 = x_PT_2
#x_final = 1000
t_initial_1 = log(x_initial_1)
t_initial_2 = log(x_initial_2)
t_final_1 = log(x_final_1)
t_final_2 = log(x_final_2)

t_q_equi_initial = log(x_PT_1)
t_q_equi_final = log(x_PT_1*4)
t_q_equi_vec = collect(t_q_equi_initial:Delta_t:t_q_equi_final)
xvec_q_equi_vec = exp.(t_q_equi_vec)
Npoints_q_equi_extra = length(t_q_equi_vec)
EquilibriumYieldQuark_extra = zeros(Npoints_q_equi_extra)
for i = 1:Npoints_q_equi_extra
    EquilibriumYieldQuark_extra[i] = Yeq(gDM, g_effective, xvec_q_equi_vec[i])
end

#First define initial conditions:
tvec_1 = collect(t_initial_1:Delta_t:t_final_1)
xvec_1 = exp.(tvec_1)
Npoints_1 = length(tvec_1)
tvec_2 = collect(t_initial_2:Delta_t:t_final_2)
xvec_2 = exp.(tvec_2)
Npoints_2 = length(tvec_2)

EquilibriumYield_1 = zeros(Npoints_1)
for i = 1:Npoints_1
    EquilibriumYield_1[i] = Yeq(gDM, g_effective, xvec_1[i]) #For testing purposes
end
EquilibriumYield_2 = zeros(Npoints_2)
for i = 1:Npoints_2
    EquilibriumYield_2[i] = Yeq(gDM, g_effective, xvec_2[i]) #For testing purposes
end

sigma_v_averaged_1 = ones(Npoints_1) #Initialise array for interpolated annihaltion xs
for i in 1:Npoints_1
    sigma_v_averaged_1[i] = sigma_1 #For testing purposes
end
sigma_v_averaged_2 = ones(Npoints_2) #Initialise array for interpolated annihaltion xs
for i in 1:Npoints_2
    sigma_v_averaged_2[i] = sigma_2 #For testing purposes
end

Wx_1 = zeros(Npoints_1)
Yx_1 = zeros(Npoints_1)
Yx_1[1] = EquilibriumYield_1[1]
Wx_1[1] = log(Yx_1[1])

Wx_2 = zeros(Npoints_2)
Yx_2 = zeros(Npoints_2)
Yx_2[1] = EquilibriumYield_2[1]
Wx_2[1] = log(Yx_2[1])

#Solution to the Boltzmann equation for the first freeze-out
for i = 2:Npoints_1
    W_old = Wx_1[i-1]
    #Wx[i] = Newton_Raphson_step(tvec[i], W_old, 0.5*BigConstant*g_star_eff_vec[i]*Delta_t*exp(-tvec[i])*sigma_v_averaged[i], g_quark, h_eff_dof(m_quark/xvec[i]))
    Wx_1[i] = Newton_Raphson_step(xvec_1[i], W_old, 0.5*BigConstant_1*g_effective*Delta_t/xvec_1[i]*sigma_v_averaged_1[i], gDM, g_effective) #For testing purposes
end
Yx_1 = exp.(Wx_1)
#Solution to the Boltzmann equation for the first freeze-out
for i = 2:Npoints_2
    W_old = Wx_2[i-1]
    #Wx[i] = Newton_Raphson_step(tvec[i], W_old, 0.5*BigConstant*g_star_eff_vec[i]*Delta_t*exp(-tvec[i])*sigma_v_averaged[i], g_quark, h_eff_dof(m_quark/xvec[i]))
    Wx_2[i] = Newton_Raphson_step(xvec_2[i], W_old, 0.5*BigConstant_2*g_effective*Delta_t/xvec_2[i]*sigma_v_averaged_2[i], gDM, g_effective) #For testing purposes
end
Yx_2 = exp.(Wx_2)

### FOPT: Squeezeout step ###
Yx_squeezeout_1 = 1.5/pi*sqrt(15*Yx_2[Npoints_1]/(2*pi*h_eff_dof(Tcrit_1)*R_pocket^3))
Yx_squeezeout_2 = 1.5/pi*sqrt(15*Yx_2[Npoints_2]/(2*pi*h_eff_dof(Tcrit_2)*R_pocket^3))

### Include a new 'freeze-out' solution stage, now for baryons
x_initial_12 = x_PT_1
x_final_12 = 1E6*m_quark_2*Lambda_dQCD_1/m_quark_1/Lambda_dQCD_2
t_initial_12 = log(x_initial_12)
t_final_12 = log(x_final_12)

t_B_equi_initial = log(x_PT_1/10)
t_B_equi_final = log(x_PT_1)
t_B_equi_vec = collect(t_B_equi_initial:Delta_t:t_B_equi_final)
xvec_B_equi_vec = exp.(t_B_equi_vec)
Npoints_B_equi_extra = length(t_B_equi_vec)
EquilibriumYieldBaryon_extra = zeros(Npoints_B_equi_extra)
for i = 1:Npoints_B_equi_extra
    EquilibriumYieldBaryon_extra[i] = Yeq_baryon(g_Baryon, h_eff_dof(m_quark_1/xvec_B_equi_vec[i]), xvec_B_equi_vec[i], alpha_1, Ndark)
end

#First define initial conditions:
tvec_12 = collect(t_initial_12:Delta_t:t_final_12)
xvec_12 = exp.(tvec_12)
Npoints_12 = length(tvec_12)

EquilibriumYieldBaryon = zeros(Npoints_12)
for i = 1:Npoints_12
    EquilibriumYieldBaryon[i] = Yeq_baryon(g_Baryon, h_eff_dof(m_quark_1/xvec_12[i]), xvec_12[i], alpha_1, Ndark)
end

sigma_v_averaged_12 = ones(Npoints_12) #Initialise array for interpolated annihaltion xs

for i in 1:Npoints_12
    sigma_v_averaged_12[i] = sigma_12 #For testing purposes
end

Wx_12 = zeros(Npoints_12)
Yx_12 = zeros(Npoints_12)
Yx_12[1] = Yx_squeezeout_1
Wx_12[1] = log(Yx_12[1])

#Solution to the Boltzmann equation for the second "freeze-out"
for i = 2:Npoints_12
    W_old = Wx_12[i-1]
    Wx_12[i] = Newton_Raphson_step_Baryon(xvec_12[i], W_old, 0.5*BigConstant_2*g_effective*Delta_t/xvec_12[i]*sigma_v_averaged_12[i], g_Baryon, g_effective, alpha_2, Ndark) #For testing purposes
end
Yx_12 = exp.(Wx_12)

### Entropy dilution due to glueball decay ###
m_glueball = 7*Lambda_dQCD_2 #Mass of the lightest 0++ glueball
x_freeze_out = GB_freeze_out_estimate(1, m_glueball, R_max) #freeze-out of dark gluons. xft = 1 is the initial guess 

#Matter-radiation equality
Y_GB = R_max/x_freeze_out #Relic yield of dark gluons
T_MR = 4/3*m_glueball*Y_GB # Matter-radiation equality temperature, after which GB dominate the energy content
x_MR = m_quark_2/T_MR

#Glueball decay
Gamma_GB = glueball_decay(m_glueball, m_quark_2, alpha_2)

#Entropy dilution
dil_fac = (1 + 1.65*g_average(1e-5, T_MR, 10)*cbrt(T_MR^4/(Gamma_GB*Mpl))^2)^(-0.75)

Yx_dilution = dil_fac*Yx_squeezeout_2

relic_abundance_Yield_limit = Relic_abundance_limit*rho_crit/(reduced_Hubble_squared*s0*Lambda_dQCD_2) 

popfirst!(xvec_12)
popfirst!(Yx_12)

plot(Lambda_dQCD_1*xvec_1/m_quark_1, EquilibriumYield_1*m_quark_1/Lambda_dQCD_1, xaxis=:log, yaxis=:log, xlabel=L"\Lambda/T", ylabel=L"Y*m_Q/\Lambda", legend=false,
 linewidth = 3, linestyle =:dash, c=:red, xlims = (1E-4, 1E7), ylims = (1E-22, 30), minorgrid = false, xticks=5, yticks = 5, label = false)
plot!(Lambda_dQCD_2*xvec_2/m_quark_2, EquilibriumYield_2*m_quark_2/Lambda_dQCD_2, linewidth = 3, linestyle =:dash, c=:blue)
plot!(Lambda_dQCD_1*xvec_q_equi_vec/m_quark_1, EquilibriumYieldQuark_extra*m_quark_1/Lambda_dQCD_1, linestyle=:dot, linewidth = 3, c=:red)
plot!(Lambda_dQCD_1*xvec_1/m_quark_1, Yx_1*m_quark_1/Lambda_dQCD_1, linewidth=3, c=:red)
plot!(Lambda_dQCD_1*xvec_12/m_quark_1, Yx_12*m_quark_1/Lambda_dQCD_1, linewidth=3, c=:red) 
plot!(Lambda_dQCD_1*xvec_B_equi_vec/m_quark_1, EquilibriumYieldBaryon_extra*m_quark_1/Lambda_dQCD_1, linewidth=3, c=:red, linestyle=:dot)
plot!(Lambda_dQCD_2*xvec_2/m_quark_2, Yx_2*m_quark_2/Lambda_dQCD_2, linewidth=3, c=:blue)
hline!([relic_abundance_Yield_limit], linestyle=:dash, c=:green, linewidth=2)
plot!([Lambda_dQCD_2*x_PT_2/m_quark_2, 1E4], [Yx_squeezeout_2*m_quark_2, Yx_squeezeout_2*m_quark_2], linewidth=3, c=:blue)
vline!([1.2], linestyle=:dashdotdot, c=:black, linewidth=5)
plot!([1E4, 1E10], [Yx_dilution*m_quark_2, Yx_dilution*m_quark_2], linewidth=3, c=:blue)
annotate!(5E-3, 1E-20, text(L"Y^Q_{eq}", :blue, :right, 12))
annotate!(5, 1E-5, text(L"Y_{eq}^Q", :red, :left, 12))
annotate!(0.5, 1E-3, text(L"Y_{eq}^B", :red, :right, 12))
annotate!(1E7, 5, text(L"\Lambda = 50"*" GeV", :black, :right, 12))
annotate!(1E7, 1E-2, text(L"m_Q/\Lambda = 5;"*L" T_c > T_{f.o.}", :red, :right, 12))
annotate!(1E7, 1E-4, text(L"m_Q/\Lambda = 10^4;"*L" T_c < T_{f.o.}", :blue, :right, 12))
annotate!(1E7, 1E-10, text(L"\Omega h^2 = 0.12", :green, :right, 12))
annotate!(1E7, 5E-19, text(L"Y_\infty(m_Q = 500"*" TeV)", :blue, :right, 10))
annotate!(1E7, 5E-21, text(L"Y_\infty(m_Q = 250"*" GeV)", :red, :right, 10))
savefig("Test.png")