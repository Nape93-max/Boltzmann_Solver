#using Plots #Plots not needed in this version of the code. Still kept for developer's purposes
using LaTeXStrings
#using LogarithmicNumbers #Not really needed

include("../FreezeOut.jl") #Include important functions 
include("../SqueezeOut.jl") #Include important functions regarding seeze-out and GB decay

function cross_section_contributions(m_N, m_L, Lambda_dQCD, ydark)
    Alpha_NN = running_coupling_from_pole(2*m_N, Lambda_dQCD, (11*Ndark-2)/3) #Alpha_DM at the scale of NN annihilation. The N flavour is active. 
    Alpha_NL = running_coupling_from_pole(m_N + m_L, Lambda_dQCD, (11*Ndark-2*2-2)/3) #Alpha_DM at the scale of NL annihilation. Both flavours are active
    Alpha_LL = running_coupling_from_pole(2*m_L, Lambda_dQCD, (11*Ndark - 2*2-2)/3) #Alpha_DM at the scale of LL annihilation
    aw_NL = running_coupling_from_scale(m_N + m_L, Mtop, alpha_W_Mtop, (19 - 4*Ndark)/6) #beta0 = 19/6 for the weak interaction below m_L
    aw_LL = running_coupling_from_scale(2*m_L, Mtop, alpha_W_Mtop, (19 - 4*Ndark)/6) #beta0 = 19/6 for the weak interaction below m_L
    aY_NL = running_coupling_from_scale(m_N + m_L, Mtop, alpha_Y_Mtop, -(41 + Ndark*4*2)/6) #beta0 = -41/6 for the U(1) interaction below m_L
    aY_LL = running_coupling_from_scale(2*m_L, Mtop, alpha_Y_Mtop, -(41 + Ndark*4*2)/6) #beta0 = -41/6 for the U(1) interaction below m_L

    aux1NL = 3*aw_NL + aY_NL # Auxiliary quantities which are defined to shorten the expressions for the cross sections a bit
    aux1LL = 3*aw_LL + aY_LL
    aux2 = 101*aw_LL*aw_LL + 8*aw_LL*aY_LL + 43*aY_LL*aY_LL
    aux3 = NdarkAdjoint*(NdarkAdjoint-1)

    sigmaNN = aux3/(8*Ndark*(NdarkAdjoint + 1))*pert_acs(Alpha_NN, m_N) #Derivation of these expression in the Mathematica notebook
    sigmaNL = (4*NdarkAdjoint*Alpha_NL*m_N*m_L + 2*NdarkAdjoint*Alpha_NL*m_L*m_L + m_N*m_N*(2*NdarkAdjoint*Alpha_NL + Ndark*aux1NL))/(32*pi)*pert_acs(ydark, Ndark*m_N*m_L)
    sigmaLL = 0.5*(2*pi*m_N*m_N*m_L*m_L*(Ndark*Ndark*(3*pi*aux2 - ydark*ydark*aux1LL) + 16*pi*Ndark*NdarkAdjoint*Alpha_LL*aux1LL + 8*pi*aux3*Alpha_LL*Alpha_LL) + m_L^4*(-2*pi*(NdarkAdjoint + 1)*ydark*ydark*aux1LL + 4*(NdarkAdjoint + 1)*ydark^4 + pi*pi*(16*NdarkAdjoint*Ndark*Alpha_LL*aux1LL + 8*aux3*Alpha_LL*Alpha_LL + 3*aux2*(NdarkAdjoint + 1))) + pi*pi*m_N^4*(16*Ndark*NdarkAdjoint*Alpha_LL*aux1LL + 8*aux3*Alpha_LL*Alpha_LL + 3*(NdarkAdjoint + 1)*aux2))*pert_acs(1, 8*pi*m_L*(m_N*m_N + m_L*m_L))/(Ndark*(NdarkAdjoint + 1))

    return sigmaNN, sigmaNL, sigmaLL
end

function cross_section(x, m_N, m_L, Lambda_dQCD, ydark, sigmaNN, sigmaNL, sigmaLL) #Calculate the averaged cross section according to Griest and Seckel with the aid of my Mathematica notebook
    EquilibriumYieldN = Yeq(g_N, h_eff_dof(m_L/m_N*x), m_L/m_N*x)
    EquilibriumYieldL = Yeq(g_L, h_eff_dof(m_L/m_N*x), m_L/m_N*x)
    EquilibriumYieldDM = EquilibriumYieldN + EquilibriumYieldL

    #=
    #Here the thermally averaged cross section is read in. 
    sigma_v_file = DataFrame(CSV.File("Sigma_eff.txt"))
    sigma_v_x_values = sigma_v_file[!,1]
    #sigma_v_y_values = sigma_v_file[!,2] #Not yet ready
    sigma_v_y_values = ones(500) #for evaluation w.o. SE and BSF

    sigma_v_averaged = ones(Npoints) #Initialise array for interpolated annihaltion xs
    sigma_v_averaged_coeffs = sigma_v_interpolation(sigma_eff, sigma_v_x_values, sigma_v_y_values) #Cubic spline fit coefficient vector (beta, gamma, delta)
    =#

    #Calculate the effective cross section according to Griest and Seckel 
    Yeq_ratios = 0.
    if EquilibriumYieldN > 0
        Yeq_ratios = EquilibriumYieldL/EquilibriumYieldN
    end
    r_GS_N = 1/(1 + Yeq_ratios) #ratios as defined in Griest and Seckel
    r_GS_L = Yeq_ratios/(1 + Yeq_ratios)

    sigma_v_averaged = (sigmaNN*r_GS_N*r_GS_N + 2*sigmaNL*r_GS_N*r_GS_N + sigmaLL*r_GS_N*r_GS_N) #Effective xs according to Griest and Seckel
end

function coannihilation_freeze_out(x_PT, m_N, m_L, Lambda_dQCD, ydark) # function that calculates the yield of the quark freeze-out
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
    BigConstant = bc_constant(m_N)

    #Define the initial conditions
    EquilibriumYieldN_initial = Yeq(g_N, h_eff_dof(m_L/m_N*xvec[1]), m_L/m_N*xvec[1])
    EquilibriumYieldL_initial = Yeq(g_L, h_eff_dof(m_L/m_N*xvec[1]), m_L/m_N*xvec[1])
    EquilibriumYieldDM_initial = EquilibriumYieldN_initial + EquilibriumYieldL_initial
    Wx = zeros(Npoints)
    Yx = zeros(Npoints)
    Yx[1] = EquilibriumYieldDM_initial
    Wx[1] = log(Yx[1])

    sigma_v_averaged = zeros(Npoints) #Calculation of the averaged cross section
    (sigma_NN, sigma_NL, sigma_LL) = cross_section_contributions(m_N, m_L, Lambda_dQCD, ydark)
    for i in 1:Npoints
        sigma_v_averaged[i] = cross_section(xvec[i], m_N, m_L, Lambda_dQCD, ydark, sigma_NN, sigma_NL, sigma_LL)
    end

    #Solution to the Boltzmann equation for the first freeze-out
    for i = 2:Npoints
        W_old = Wx[i-1]
        Wx[i] = Newton_Raphson_step(tvec[i], W_old, 0.5*BigConstant*Delta_t*sigma_v_averaged[i], g_N, h_eff_dof(m_N/xvec[i]))
    end
    Yx = exp.(Wx)
    return Yx[Npoints-1]
end

function glueball_decay(m_glueball, Lambda_dQCD, m_N, m_L, ydark)
    Alpha_DM_GB_decay = running_coupling_from_pole(m_glueball, Lambda_dQCD, 11*Ndark/3) #Dark gauge coupling at the mass scale of the GBs
    Alpha_weak_GB_decay = running_coupling_from_scale(m_glueball, Mtop, alpha_W_Mtop, 19/6) #Weak gauge coupling at the mass scale of the GBs
    Alpha_Y_GB_decay = running_coupling_from_scale(m_glueball, Mtop, alpha_Y_Mtop, -41/6) #Hypercharge gauge coupling at the mass scale of the GBs
    decay_const_GB = 3.06*m_glueball^3/(4*pi*Alpha_DM_GB_decay) #decay constant of the gluon after Juknevich. 
    Gamma_GB_dim8 = Alpha_DM_GB_decay*Alpha_DM_GB_decay/(8*pi*m_L^8)*1/3600*m_glueball^3*(decay_const_GB)^2*(Alpha_Y_GB_decay*Alpha_Y_GB_decay/4 + 1.5*Alpha_weak_GB_decay*Alpha_weak_GB_decay) #dim 8 Glueball decay rate after Juknevich.
    Gamma_GB_dim6 = (ydark*ydark*Alpha_DM_GB_decay*decay_const_GB/(m_L*m_N))^2/(72*pi^3*m_glueball)  #dim 6 Glueball decay rate after Juknevich.
    Gamma_GB = Gamma_GB_dim6 + Gamma_GB_dim8
end

#Define constant physics parameters
const Ndark = 3. #Dark SU(N) parameter N
const NdarkAdjoint = Ndark*Ndark-1

#Running of SU(2)L gauge coupling
const MZ = 91.1876 # Z boson mass in GeV
const Mtop = 173 # Top quark mass in GeV
const alpha_W_MZ = 0.0342556 # Source?!
const alpha_W_Mtop = 0.0334 # The weak gauge coupling at the top quark mass (All from arxiv: 1307.3536)
const alpha_Y_Mtop = 0.0102 # Source: arxiv: 1307.3536

#Constant parameters for entropy dilution
const R_max = 2.5E-4 #highest possible entropy ratio after the PT
const BBN_lifetime = 6.58*1E-25 #Lower bound on glueball decay rate.

num_scales = 100
num_masses = 100
num_deltas = 10
num_Yukawas = 10
num_parameter_points = num_scales*num_masses*num_deltas*num_Yukawas
array_scales = 10.0.^collect(range(0, 7, num_scales)) #10.0.^collect(range(0, 7, length = num_scales)) 
array_masses = 10.0.^collect(range(2, 4, num_masses)) 
array_deltas = 10.0.^collect(range(0, 2, num_deltas)) 
array_Yukawas = 10.0.^collect(range(-7, -1, num_Yukawas)) 

const g_N = 4*Ndark #degeneracy of the Dirac quark N: (Spin x Particle-Antiparticle) x DarkColour
const g_L = 4*Ndark*2 #degeneracy of the Dirac quark L: (Spin x Particle-Antiparticle) x DarkColour x weak multiplicity

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

Threads.@threads for (i,j,k,l) in collect(Iterators.product(1:num_scales, 1:num_masses, 1:num_deltas, 1:num_Yukawas)) 
    big_ind = Int(i + num_scales*((j-1) + num_masses*((k-1) + num_deltas*(l-1))))

    Lambda_dQCD = array_scales[i]
    m_N = array_masses[j]*Lambda_dQCD
    mass_delta = array_deltas[k]
    m_L = (mass_delta+1)*m_N 
    ydark = array_Yukawas[l] #dark Yukawa coupling. No running implemented. 

    Alpha_DM = running_coupling_from_pole(m_N, Lambda_dQCD, 11*Ndark/3)
    AlphaDM_data_vec[big_ind] = Alpha_DM

    Tcrit = 0.63*Lambda_dQCD #Temperature of the phase transition
    x_PT = m_N/Tcrit 
    xPT_data_vec[big_ind] = x_PT

    #Constants of FOPT (w/o factors of 1/Lambda for numerical convenience, cancel in the squeezeout step!)
    R_pocket = pocket_radius(Lambda_dQCD)
    RPocket_data_vec[big_ind] = R_pocket

    Yfo = coannihilation_freeze_out(x_PT, m_N, m_L, Lambda_dQCD, ydark)
    Yfo_data_vec[big_ind] = Yfo

    ### FOPT: Squeezeout step ###
    Yx_squeezeout = 1.5/pi*sqrt(15*Yfo/(2*pi*h_eff_dof(Tcrit)*R_pocket^3))
    Ysqo_data_vec[big_ind] = Yx_squeezeout

    ### Entropy dilution due to glueball decay ###
    m_glueball = 7*Lambda_dQCD #Mass of the lightest 0++ glueball
    x_freeze_out = GB_freeze_out_estimate(1, m_glueball, R_max) #freeze-out of dark gluons
    xGBfo_data_vec[big_ind] = x_freeze_out

    #Matter-radiation equality
    Y_GB = R_max/x_freeze_out #Relic yield of dark gluons
    T_MR = 4/3*m_glueball*Y_GB # Matter-radiation equality temperature, after which GB dominate the energy content
    TMReq_data_vec[big_ind] = T_MR
    x_MR = m_N/T_MR

    #GB decay
    Gamma_GB = glueball_decay(m_glueball, Lambda_dQCD, m_N, m_L, ydark)
    Gamma_GB_data_vec[big_ind] = Gamma_GB

    #Entropy dilution
    dil_fac = (1 + 1.65*g_average(1e-5, T_MR, 10)*cbrt(T_MR^4/(Gamma_GB*Mpl))^2)^(-0.75)
    #x_dilution = x_PT*100
    Yx_dilution = dil_fac*Yx_squeezeout
    dil_fac_data_vec[big_ind] = dil_fac

    ### Calculation of the relic abundance ###
    x_today = m_N/T0
    Omega_relic = reduced_Hubble_squared*Yx_dilution*s0*m_N/rho_crit
    Omegah2_data_vec[big_ind] = Omega_relic
end

#Initialise final output file with data
results_file = open("NplusL_model_scan.csv", "w")
IOStream("NplusL_model_scan.csv")
write(results_file, "mN/GeV, Lambda/GeV, x_PT, Alpha(mN), delta, mL/GeV, ydark, RPocket/Lambda, Yfo, Ysqo, xGBfo, TMReq/GeV, Gamma_GB/GeV, dilution_factor, Omegah2\n")

#Write out the final data
for i in 1:num_scales
    for j in 1:num_masses
        for k in 1:num_deltas
            for l in 1:num_Yukawas
                local big_ind = Int(i + num_scales*((j-1) + num_masses*((k-1) + num_deltas*(l-1))))
                write(results_file, join((array_masses[j]*array_scales[i], array_scales[i], xPT_data_vec[big_ind], AlphaDM_data_vec[big_ind], array_deltas[k], (1 + array_masses[j]*array_scales[i])*array_deltas[k], array_Yukawas[l], RPocket_data_vec[big_ind], Yfo_data_vec[big_ind], Ysqo_data_vec[big_ind], xGBfo_data_vec[big_ind], TMReq_data_vec[big_ind], Gamma_GB_data_vec[big_ind], dil_fac_data_vec[big_ind], Omegah2_data_vec[big_ind]),","),"\n")
            end
        end
    end
end

close(results_file)