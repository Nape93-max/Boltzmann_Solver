include("NplusL_Header.jl") #Include important functions 

num_scales = 100
num_masses = 100
num_deltas = 10
num_Yukawas = 10
num_parameter_points = num_scales*num_masses*num_deltas*num_Yukawas
array_scales = 10.0.^collect(range(-3, 7, num_scales)) #10.0.^collect(range(0, 7, length = num_scales)) 
array_masses = 10.0.^collect(range(1.2, 4, num_masses)) 
array_deltas = 10.0.^collect(range(0, 2, num_deltas)) 
array_Yukawas = 10.0.^collect(range(-7, 0, num_Yukawas)) 

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

    Tcrit = 1.91*Lambda_dQCD #Temperature of the phase transition
    x_PT = m_N/Tcrit 
    xPT_data_vec[big_ind] = x_PT

    #Constants of FOPT (w/o factors of 1/Lambda for numerical convenience, cancel in the squeezeout step!)
    R_pocket = pocket_radius(Lambda_dQCD)
    RPocket_data_vec[big_ind] = R_pocket

    #FreezeOut
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
    x_MR = m_N/T_MR
    TMReq_data_vec[big_ind] = T_MR

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