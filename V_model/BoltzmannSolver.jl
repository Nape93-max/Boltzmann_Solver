#using Plots
using LaTeXStrings
#using LogarithmicNumbers #Not really needed

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

#Running of SU(2)L gauge coupling
const MZ = 91.1876 # Z boson mass in GeV
const Mtop = 173 # Top quark mass in GeV
const alpha_W_MZ = 0.0342556 # Source?!
const alpha_W_Mtop = 0.0334 # The weak gauge coupling at the top quark mass (All from arxiv: 1307.3536)

#Constant parameters for entropy dilution
const R_max = 2.5E-4 #highest possible entropy ratio after the PT
const BBN_lifetime = 6.58*1E-25 #Lower bound on glueball decay rate.

#Details on the parameter scan
num_scales = 100
num_masses = 100
num_parameter_points = num_scales*num_masses
array_scales = 10.0.^collect(range(0, 7, length = num_scales)) #0 - 7 
array_masses = 10.0.^collect(range(2, 4, length = num_masses)) # 2 - 4

const g_quark = 4*Ndark*3 #degeneracy of the Dirac quark: (Spin x Particle-Antiparticle) x DarkColour x weak multiplicity

#initialise arrays of quantities that will be calculated and later written into results file
xPT_data_vec = zeros(num_parameter_points)
AlphaDM_data_vec = zeros(num_parameter_points)
RPocket_data_vec = zeros(num_parameter_points)
Yfo_data_vec = zeros(num_parameter_points)
Ysqo_data_vec = zeros(num_parameter_points)
xGBfo_data_vec = zeros(num_parameter_points)
TMReq_data_vec = zeros(num_parameter_points)
Gamma_GB_data_vec = zeros(num_parameter_points)
dil_fac_data_vec = zeros(num_parameter_points)
Omegah2_data_vec = zeros(num_parameter_points)

Threads.@threads for (i,j) in collect(Iterators.product(1:length(array_scales), 1:length(array_masses))) #Masterloop of parameter scans
   big_ind = Int(i + (j-1)*num_scales)

    #Define physics parameters of the model (parameters in the loop)
    Lambda_dQCD = array_scales[i]
    m_quark = array_masses[j]*Lambda_dQCD
    Alpha_DM = running_coupling_from_pole(2*m_quark, Lambda_dQCD, 11*Ndark/3-2*3/3) #At the annihilation scale, the quark is active
    Alpha_dark = running_coupling_from_pole(m_quark, Lambda_dQCD, 11*Ndark/3-2*3/3) #At the annihilation scale, the quark is active #dark gauge coupling at the mass scale m_quark
    AlphaDM_data_vec[big_ind] = Alpha_dark

    BigConstant = bc_constant(m_quark)
    sigma0 = cross_section(m_quark, Alpha_DM)
    #Lambda_dQCD = Landau_pole(m_quark, Alpha_DM, 11) #beta0 = 11*Nc/3
    Tcrit = 0.63*Lambda_dQCD #Temperature of the phase transition
    x_PT = m_quark/Tcrit 
    xPT_data_vec[big_ind] = x_PT

    #Constants of FOPT (w/o factors of 1/Lambda for numerical convenience, cancel in the squeezeout step!)
    R_pocket = pocket_radius(Lambda_dQCD)
    RPocket_data_vec[big_ind] = R_pocket
    
    #Quark freeze-out
    Yfo = quark_freeze_out(x_PT, m_quark, sigma0, BigConstant, g_quark)
    Yfo_data_vec[big_ind] = Yfo

    ### FOPT: Squeezeout step ###
    Yx_squeezeout = 1.5/pi*sqrt(15*Yfo/(2*pi*h_eff_dof(Tcrit)*R_pocket^3))
    Ysqo_data_vec[big_ind] = Yx_squeezeout

    ### Entropy dilution due to glueball decay ###
    m_glueball = 7*Lambda_dQCD #Mass of the lightest 0++ glueball
    x_freeze_out = GB_freeze_out_estimate(1, m_glueball, R_max) #freeze-out of dark gluons. xft = 1 is the initial guess 
    xGBfo_data_vec[big_ind] = x_freeze_out

    #Matter-radiation equality
    Y_GB = R_max/x_freeze_out #Relic yield of dark gluons
    T_MR = 4/3*m_glueball*Y_GB # Matter-radiation equality temperature, after which GB dominate the energy content
    TMReq_data_vec[big_ind] = T_MR
    x_MR = m_quark/T_MR

    #Glueball decay
    Gamma_GB = glueball_decay(m_glueball, m_quark, Alpha_DM)
    Gamma_GB_data_vec[big_ind] = Gamma_GB

    #Entropy dilution
    dil_fac = (1 + 1.65*g_average(1e-5, T_MR, 10)*cbrt(T_MR^4/(Gamma_GB*Mpl))^2)^(-0.75)
    Yx_dilution = dil_fac*Yx_squeezeout
    dil_fac_data_vec[big_ind] = dil_fac

    ### Calculation of the relic abundance ###
    x_today = m_quark/T0
    Omega_relic = reduced_Hubble_squared*Yx_dilution*s0*m_quark/rho_crit
    Omegah2_data_vec[big_ind] = Omega_relic
end

#Initialise final output file with data
results_file = open("V_model_scan.csv", "w")
IOStream("V_model_scan.csv")
write(results_file, "m/GeV, Lambda/GeV, x_PT, Alpha(m), RPocket/Lambda, Yfo, Ysqo, xGBfo, TMReq/GeV, Gamma_GB/GeV, dilution_factor, Omegah2\n")
for i in 1:num_scales #Write out data
    for j in 1:num_masses
        local big_ind = Int(i + (j-1)*num_scales)
        write(results_file, join((array_masses[j]*array_scales[i], array_scales[i], xPT_data_vec[big_ind], AlphaDM_data_vec[big_ind], RPocket_data_vec[big_ind], Yfo_data_vec[big_ind], Ysqo_data_vec[big_ind], xGBfo_data_vec[big_ind], TMReq_data_vec[big_ind], Gamma_GB_data_vec[big_ind], dil_fac_data_vec[big_ind], Omegah2_data_vec[big_ind]),","),"\n")
    end
end
close(results_file)

#=
### Here the plotting business starts ###

ytics = 10.0.^collect(range(-30, -1, step=1))
xtics = 10.0.^collect(range(log10(x_initial), log10(x_today), step=1))

x_vec_squeezeout = range(x_PT, length=100, stop=x_dilution)
y_vec_squeezeout = Yx_squeezeout*ones(100)
x_vec_dilution = range(x_dilution, length=100, stop=x_today)
y_vec_dilution = Yx_dilution*ones(100)

#plot(sigma_v_x_values, sigma_v_y_values, xaxis=:log, yaxis=:log)
#savefig("sigma_v_data.png")

plot(xvec, [EquilibriumYield, Yx], title="WIMP freeze-out", label=[L"Y_{eq}(x)" L"Y(x)"], yticks = ytics, xticks = xtics, minorticks = 10, minorgrid = true, xlabel="x = m/T", ylabel="Y(x)", xaxis=:log, yaxis=:log, xlims = (x_initial, x_final*1000), ylims = (1E-30, 1E-1))
plot!(x_vec_squeezeout, y_vec_squeezeout, label = L"Y_S")
plot!(x_vec_dilution, y_vec_dilution, label = L"Y_\infty")
plot!([x_MR], seriestype = :vline, label = L"x_{MR}")
plot!([x_PT], seriestype = :vline, label = L"x_\Lambda")
plot!([x_dilution], seriestype = :vline, label = L"x_{GB}")
savefig("FreezeOut.png")

plot(xvec, Yx, title="Relic Yield", minorticks = 10, minorgrid = true, xlabel=L"x = m/T", ylabel=L"Y(x)", label = L"Y(x)", yticks = ytics, xaxis=:log, yaxis=:log, xticks = xtics, xlims = (x_initial, x_final), ylims = (1E-16, 1E-1))
savefig("Y_plot.png")
=#