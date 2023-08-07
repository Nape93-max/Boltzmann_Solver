using Plots #Plots not needed in this version of the code. Still kept for developer's purposes
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

function cross_section(x, m_N, m_L, sigmaNN, sigmaNL, sigmaLL) #Calculate the averaged cross section according to Griest and Seckel
    # with the aid of my Mathematica notebook
    EquilibriumYieldN = Yeq(g_N, h_eff_dof(m_N/x), x)
    EquilibriumYieldL = Yeq(g_L, h_eff_dof(m_N/x), m_L/m_N*x)
    #EquilibriumYieldDM = EquilibriumYieldN + EquilibriumYieldL

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
    #r_GS_L = Yeq_ratios/(1 + Yeq_ratios)

    return sigmaNN*r_GS_N*r_GS_N + 2*sigmaNL*r_GS_N*r_GS_N + sigmaLL*r_GS_N*r_GS_N #Effective xs according to Griest and Seckel
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

    mass_ratio = m_L/m_N
    temperature_vec = m_N./xvec

    g_star_eff_vec = eff_dof_sqrt.(temperature_vec) #Effective degrees of freedom
    h_eff_dof_vec = h_eff_dof.(temperature_vec) #Effective entropic dofs
    BigConstant = bc_constant(m_N)

    #Define the initial conditions
    EquilibriumYieldDM_initial = Yeq_DM_coannihilation(g_N, g_L, h_eff_dof_vec[1], mass_ratio, xvec[1])
    Wx = zeros(Npoints)
    Yx = zeros(Npoints)
    Yx[1] = EquilibriumYieldDM_initial
    Wx[1] = log(Yx[1])

    mod_BC = 0.5*BigConstant*Delta_t #modified Big Constant
    sigma_v_averaged = zeros(Npoints) #Calculation of the averaged cross section
    (sigma_NN, sigma_NL, sigma_LL) = cross_section_contributions(m_N, m_L, Lambda_dQCD, ydark)

    #=  #Interpolation has proven to be slower than direct computation!
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
    #### END interpolation of sigma_v
    =#

    for i in 1:Npoints
        #sigma_v_averaged[i] = sigma_v_cspline(xvec[i], x_raw_data_vec, sigma_v_raw_data, sigma_v_raw_data_first_entry, interpolation_coeffs)
        sigma_v_averaged[i] = cross_section(xvec[i], m_N, m_L, sigma_NN, sigma_NL, sigma_LL)
    end

    #Solution to the Boltzmann equation for the first freeze-out
    for i = 2:Npoints
        W_old = Wx[i-1]
        Wx[i] = Newton_Raphson_step_coannihilation(xvec[i], W_old, mod_BC*g_star_eff_vec[i]*sigma_v_averaged[i]/xvec[i], g_N, g_L, h_eff_dof_vec[i], mass_ratio)
    end
    Yx = exp.(Wx)
    return Yx[Npoints-1]
end

function Newton_Raphson_step_coannihilation(x, W_old, cs, g_deg_1, g_deg_2, gS, mass_ratio) #Method to calculate the next W=log(Y) point in the implicit backward Euler method.
    # x is the x_(n+1) time step, W_old is W_n and cs is Delta_t*lambda(t_(n+1))*geff(t_(n+1)). g_deg is the degeneracy of the annihilating particles.
    # gS are the effective entropy dofs. This is a particular function for the coannihialtion scenario
    W_try_initial = W_old
    W_new = Newton_Raphson_iteration_coannihilation(x, W_old, cs, W_try_initial, g_deg_1, g_deg_2, gS, mass_ratio)
    diff = abs(log(abs(W_new/W_try_initial)))
    while diff > 1E-2 
        W_try = deepcopy(W_new);
        W_new = Newton_Raphson_iteration_coannihilation(x, W_old, cs, W_try, g_deg_1, g_deg_2, gS, mass_ratio)
        diff = abs(log(abs(W_new/W_try)))
    end
    if isnan(W_new)
        println("ALARM: W_new in function Newton_Raphson is NaN. Maybe relax bound on diff in the while loop as a fix.")
    end
    return W_new
end

function Newton_Raphson_iteration_coannihilation(x, W_old, cs, W_previous, g_deg_1, g_deg_2, gS, mass_ratio) #Does one NR-step to calculate a trial W^(i+1)_(n+1).
    # in the coannihilation scenario. W_previous = W^(i)_(n+1)
    A = exp(W_previous);
    B = Yeq_DM_coannihilation(g_deg_1, g_deg_2, gS, mass_ratio, x)^2/A 
    W_next = W_previous - (W_previous - W_old + cs*(A - B))/(1 + cs*(A + B))
    return W_next
end

function Yeq_DM_coannihilation(g_deg_1, g_deg_2, gS, mass_ratio, x) #Calculates the total N.R. DM equilibrium yield
    # in a coannihilation scenario with 2 species with degeneracies g_deg_1 and g_deg_2 with gS entropic dofs.
    # mass_ratio is m_L/m_N and x is m_N/T 
    EquilibriumYieldN = Yeq(g_deg_1, gS, x)
    EquilibriumYieldL = Yeq(g_deg_2, gS, mass_ratio*x)
    return EquilibriumYieldN + EquilibriumYieldL
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
const Relic_abundance_limit = 0.12 

const g_N = 4*Ndark #degeneracy of the Dirac quark N: (Spin x Particle-Antiparticle) x DarkColour
const g_L = 4*Ndark*2 #degeneracy of the Dirac quark L: (Spin x Particle-Antiparticle) x DarkColour x weak multiplicity