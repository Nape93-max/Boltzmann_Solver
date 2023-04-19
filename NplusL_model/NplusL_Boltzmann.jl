#using Plots #Plots not needed in this version of the code. Still kept for developer's purposes
using LaTeXStrings
using SpecialFunctions
using CSV, DataFrames
using QuadGK
#using LogarithmicNumbers #Not really needed

function bc_constant(m)
    c = sqrt(pi/45)*Mpl*m
end

function Newton_Raphson_step(t, W_old, cs, g_deg, gS) #Method to calculate the next W=log(Y) point in the implicit backward Euler method.
    # t = log(x) is the t_(n+1) time step, W_old is W_n and cs is Delta_t*lambda(t_(n+1))*geff(t_(n+1)). g_deg is the degeneracy of the annihilating particles.
    # gS are the effective entropy dofs. 
    W_try_initial = W_old
    W_new = Newton_Raphson_iteration(t, W_old, cs, W_try_initial, g_deg, gS)
    diff = abs(log(abs(W_new/W_try_initial))) #before, with less numerical precision: abs(log(W_new/W_try_initial))
    while diff > 1E-2
        W_try = copy(W_new);
        W_new = Newton_Raphson_iteration(t, W_old, cs, W_try, g_deg, gS)
        diff = abs(log(abs(W_new/W_try)))
    end
    if isnan(W_new)
        println("ALARM: W_new in function Newton_Raphson is NaN. Maybe relax bound on diff in the while loop as a fix.")
    end
    return W_new
end

function Newton_Raphson_iteration(t, W_old, cs, W_previous, g_deg, gS) #Does one NR-step to calculate a trial W^(i+1)_(n+1). W_previous = W^(i)_n
    A = exp(W_previous);
    B = Yeq(g_deg, gS, exp(t))^2/A
    W_next = W_previous - (W_old - W_previous + cs*(A - B))/(1 + cs*(A + B))
end

function pert_acs(alpha, m)
    h = alpha/m
    sigma = pi*h*h
end

function Yeq(g_deg, gS, x) #Calculates the equilibrium yield of species with degeneracy g_deg at time x with gS relativistic entropic d.o.f.
    if x < 100
        11.25*g_deg*(x/(pi*pi))^2/gS*besselk(2,x)
    else
        90/(2*pi)^3.5*g_deg/gS*x*sqrt(x)*exp(-x) #This is the non-relativistic version of the line above
    end
end

function eff_dof_sqrt(T)
    if T <= first(temperatures)
        g = first(g_star_eff_sqrt)
    elseif T >= last(temperatures)
        g = last(g_star_eff_sqrt)
    else
        ind = findfirst(temperatures -> temperatures > T, temperatures);
        To = temperatures[ind-1]
        Tn = temperatures[ind]
        go = g_star_eff_sqrt[ind-1]
        gn = g_star_eff_sqrt[ind]
        g = go + (gn - go)*(T-To)/(Tn-To) # simple linear interpolation
    end
end

function h_eff_dof(T)
    if T <= first(temperatures)
        h = first(heff)
    elseif T >= last(temperatures)
        h = last(heff)
    else
        ind = findfirst(temperatures -> temperatures > T, temperatures);
        To = temperatures[ind-1]
        Tn = temperatures[ind]
        ho = heff[ind-1]
        hn = heff[ind]
        h = ho + (hn - ho)*(T-To)/(Tn-To) # simple linear interpolation
    end
end

function g_eff_dof(T)
    if T <= first(temperatures)
        h = first(geff)
    elseif T >= last(temperatures)
        h = last(geff)
    else
        ind = findfirst(temperatures -> temperatures > T, temperatures);
        To = temperatures[ind-1]
        Tn = temperatures[ind]
        ho = geff[ind-1]
        hn = geff[ind]
        h = ho + (hn - ho)*(T-To)/(Tn-To) # simple linear interpolation
    end
end

function Landau_pole(mu, alpha, beta0) #This calculates the Landau pole for a theory with coefficient beta0 given a coupling alpha(mu) at another scale mu. 
    Lambda = mu*exp(-2*pi/(beta0*alpha))
end

function running_coupling_from_pole(Q, Lambda, beta0) #Running coupling at the scale Q given a Landau pole Lambda with a theory with a beta0 coefficient
    return 2*pi/(beta0*log(Q/Lambda))
end

function running_coupling_from_scale(Q, mu, alphamu, beta0) #Running coupling at the scale Q given the running alphamu at another scale mu with a beta0 coefficient. 
    return alphamu/(1 + alphamu*beta0*log(Q/mu)/(2*pi))
end

function sigma_v_cspline(x, x_input, y_input, cxs, coeffs) #Cubic spline interpolation step with the previously calculated coefficient vectors beta, gamma and delta encoded in coeffs. x is the argument, x(y)_input is the input x(y)-data and cxs is the constant part of the cross section
    n = length(x_input)
    xvals = log.(x_input) #logarithmic interpolation
    yvals = log.(y_input)
    if x > last(x_input)
        ind = n
    else
        ind = findfirst(x_input -> x_input > x, x_input) #find relevant interval
    end

    beta = coeffs[ind] #coefficients beta_i, gamma_i and delta_i
    gamma = coeffs[ind + n]
    delta = coeffs[ind + 2*n]

    arg = log(x) - xvals[ind] # x-x_i
    p = cxs*exp(yvals[ind] + arg*(beta + arg*(gamma + arg*delta))) #Cubic spline step
end

function Gauss_elim_spline(n, mu, lambda, d) #Calculates the solution vector M used in the function sigma_v_interpolation via Gaussian elimination of a tridiagonal matrix. n is the dimension of the matrix, mu and lambda are arrays of numbers appearing in the matrix (calculated in the mother function). d is the array of the inhomogeneity. 
    #Initialisation of quantities
    M = zeros(n)
    M0 = zeros(n)
    alpha = 1
    q = zeros(n)
    beta = zeros(n-1)
    rho = zeros(n-1)
    rho_q = zeros(n-1)

    #Calculation of auxiliary vectors beta, rho and rho_q
    beta[1] = 2 - 0.5*lambda[1]*mu[2]
    rho[1] = d[2] - 0.5*d[1]*mu[2]
    rho[1] = -0.5*alpha*mu[2] #rho_q_2 = u(2) - u(1)*mu(2)/2
    for i = 2:n-1
        beta[i] = 2 - lambda[i]*mu[i+1]/beta[i-1]
        rho[i] = d[i] - rho[i-1]*mu[i+1]/beta[i] 
        rho_q[i] = -rho_q[i-1]*mu[i+1]/beta[i]
    end
    rho_q[n-1] += lambda[n]

    #Calculation of M0 and q 
    M0[n] = rho[n-1]/beta[n-1]
    q[n] = rho_q[n-1]/beta[n-1]
    for i in reverse(3:n)
        M0[i-1] = (rho[i-2] - lambda[i-1]*M0[i])/beta[i-2]
        q[i-1] = (rho_q[i-2] - lambda[i-1]*q[i])/beta[i-2]
    end
    M0[1] = (d[1] - lambda[1]*M0[2])/beta[1]
    q[1] = (alpha - lambda[1]*q[2])/beta[1]

    #Final step
    vdotM0 = M0[1] + M0[n]*mu[1]/alpha
    vdotq = q[1] + mu[1]*q[n]/alpha
    convenient_constant = vdotM0/(1 + vdotq)
    for i = 1:n
        M[i] = M0[i] - q[i]*convenient_constant
    end
    return M
end

function sigma_v_interpolation(cs, x_input, y_input) # Interpolation of sigma_v with a cubic spline. x(y)_input are the vectors containing input data of the xs. cs is a constant part of the xs. The output is a 3*n array of (beta; gamma; delta) coefficients for the spline.
    n = length(x_input) #read in data in logarithmic scaling for improved accuracy
    xvals = log.(x_input)
    yvals = log.(y_input)

    beta = zeros(n) # initialise output vectors
    gamma = zeros(n)
    delta = zeros(n)

    A = zeros(n) # Initialise auxiliary vectors for the spline

    h = prepend!(diff(xvals), xvals[1] - xvals[n]) #Set matrix entries for the spline matrix problem
    hfirst = h[1] #auxiliary variable
    y_diff_vec = prepend!(diff(yvals), yvals[1] - yvals[n]) #Auxiliary expr appearing in the definition of the inhomogeneity d_i
    lambda = ones(n) 
    d = zeros(n)
    mu = zeros(n)

    for i = 1:n-1
        lambda[i] = 1/(1 + h[i]/h[i+1])
        mu[i] = 1 - lambda[i] 
        d[i] = 6/(h[i] + h[i+1])*(y_diff_vec[i+1]/h[i+1] - y_diff_vec[i]/h[i])
    end
    lambda[n] = 1/(1 + h[n]/hfirst)
    mu[n] = 1 - lambda[n]
    d[n] = 6/(hfirst + h[n])*(y_diff_vec[1]/hfirst - y_diff_vec[n]/h[n])

    ### Solution via Gauss elimination
    M = Gauss_elim_spline(n, mu, lambda, d)
    ###
    M_diff_vec = append!(diff(M), M[1] - M[n])

    for i = 1:n-1 #Now calculate the coefficients beta_i, gamma_i and delta_i from the solution of the matrix problem
        hind = h[i+1]
        A[i] = y_diff_vec[i+1]/hind - hind/6*M_diff_vec[i]
        beta[i] = A[i] - M[i]*hind/2
        delta[i] = M_diff_vec[i]/(6*hind)
    end
    A[n] = y_diff_vec[1]/hfirst - hfirst/6*M_diff_vec[n] #Declare the boundary terms
    beta[n] = A[n] - M[n]*hfirst/2
    gamma = 0.5.*M 
    delta[n] = M_diff_vec[n]/(6*hfirst)

    coeffs = zeros(3*n) #Initialising the final result vector
    for i = 1:n
        coeffs[i] = beta[i]
        coeffs[i + n] = gamma[i]
        coeffs[i + 2*n] = delta[i]
    end
    return coeffs
end

function aux_func_GB_decay(xfo, m_GB, R)
    sigma32 = (4*pi)^3/Ndark^6/m_GB
    lhs = log(xfo)*(5/2) + 2*xfo - log(h_eff_dof(m_GB/xfo)*R/(180*pi)*(Mpl*sigma32/sqrt(4/5*pi^3*g_eff_dof(m_GB/xfo)))^(3/2))
    return lhs
end

function GB_freeze_out_estimate(xft, m_GB, R) #Numerical estimate of the freeze-out x with the secant method. xft is the initial guess, m is the DM mass and R is the entropy ratio. acs is the constant annihilation xs. cs_xvalues are the x_values of the fully dressed acs. cs_yvalues are the corresponding acs values. coeffs are the spline interpolation coefficients.
    xf0 = xft #first try
    xf1 = xf0 - aux_func_GB_decay(xf0, m_GB, R)*2*0.001/(aux_func_GB_decay(xf0 + 0.001, m_GB, R)-aux_func_GB_decay(xf0 - 0.001, m_GB, R))
    diff = abs(xf1 - xf0)
    if diff < 1E-4
        xf2 = xf1
    else
        while diff > 1E-4
            xf2 = xf1 - aux_func_GB_decay(xf1, m_GB, R)*(xf1 - xf0)/(aux_func_GB_decay(xf1, m_GB, R) - aux_func_GB_decay(xf0, m_GB, R))
            diff = abs(xf2 - xf1)
            xf0 = copy(xf1)
            xf1 = copy(xf2)
        end
    end
    return xf2
end

function g_average(ui, Ti, uf) #Calculates the averaged eff dofs of entropy needed in the entropy dilution
    integral_denom, err_denom = quadgk(u->u^(2/3)*exp(-u), ui, uf, rtol=1e-8)
    integral_num, err_num = quadgk(u->u^(2/3)*exp(-u)*h_eff_dof(T_convert_u(u, ui, Ti))^(1/3), ui, uf, rtol=1e-8)
    return integral_num/integral_denom
end

function T_convert_u(u, ui, Ti) #Converts the variable u used in the g-average into temperatures
    T = Ti*(ui/u)^(2/3)
end

function rolled_index(big_ind, num_array) #This function converts a rolled index ind into the subindices small_inds = (n1, n2, ...). The maximally allowed values for the ni are stored in num_array
    # The formula that holds is 
    # ind = n1 + (n2-1)*L1 + (n3-1)*L1*L2 + ... + (nd-1)*L1*L2*...*L_(d-1)
    # Via basic arithmetics, this relation can be inverted for the ni. 
    d = length(num_array)
    small_inds = ones(d)
    
    new_array = deepcopy(num_array)
    new_ind = deepcopy(big_ind)
    
    h1 = prod(new_array) #It works, believe me.
    for j in d:-1:1 
        if new_ind == h1 || new_ind == 0
            for k in 1:j
                small_inds[k] = num_array[k]
            end
            break
        elseif j == 1
            small_inds[j] = new_ind
        else
            pop!(new_array)
            h1 = prod(new_array)
            h2 = ceil(new_ind/h1) 
            if new_ind > h1
               small_inds[j] = h2
               new_ind -= h1*(h2-1)
            end
        end
    end
    return Int.(small_inds)
end

#Define constant physics parameters
const Ndark = 3. #Dark SU(N) parameter N
const NdarkAdjoint = Ndark*Ndark-1
const Mpl = 1.221E19
const H0 = 1.447E-42
const T0 = 2.35E-13
const rho_crit = 3.724E-47
const s0 = 2.225E-38
const reduced_Hubble_squared = 0.67*0.67

#Read and define degrees of freedom 
dof_file = DataFrame(CSV.File("DegreesOfFreedom.txt"))
const temperatures = dof_file[!,1]
const geff = dof_file[!,2]
const heff = dof_file[!,3]
const g_star_eff_sqrt = dof_file[!,4]
const num_of_g_points = length(temperatures)

#Running of SU(2)L gauge coupling
const MZ = 91.1876 # Z boson mass in GeV
const Mtop = 173 # Top quark mass in GeV
const alpha_W_MZ = 0.0342556 # Source?!
const alpha_W_Mtop = 0.0334 # The weak gauge coupling at the top quark mass (All from arxiv: 1307.3536)
const alpha_Y_Mtop = 0.0102 # Source: arxiv: 1307.3536

#Constant parameters for entropy dilution
const R_max = 2.5E-4 #highest possible entropy ratio after the PT
const BBN_lifetime = 6.58*1E-25 #Lower bound on glueball decay rate.

num_scales = 2
num_masses = 2
num_deltas = 2
num_Yukawas = 2
num_of_parameter_points = num_scales*num_masses*num_deltas*num_Yukawas
array_scales = 10.0.^collect(range(0, 7, num_scales)) #10.0.^collect(range(0, 7, length = num_scales)) 
array_masses = 10.0.^collect(range(2, 4, num_masses)) 
array_deltas = 10.0.^collect(range(0, 2, num_deltas)) 
array_Yukawas = 10.0.^collect(range(-7, -1, num_Yukawas)) 

const g_N = 4*Ndark #degeneracy of the Dirac quark N: (Spin x Particle-Antiparticle) x DarkColour
const g_L = 4*Ndark*2 #degeneracy of the Dirac quark L: (Spin x Particle-Antiparticle) x DarkColour x weak multiplicity

#Initialise final output file with data
results_file = open("NplusL_model_scan.csv", "w")
IOStream("NplusL_model_scan.csv")
write(results_file, "mN/GeV, Lambda/GeV, mN/Lambda, Alpha(mN), delta, mL/GeV, ydark, RPocket/Lambda, Yfo, Ysqo, xGBfo, TMReq/GeV, Gamma_GB/GeV, dilution_factor, Omegah2\n")

Threads.@threads for (i,j,k,l) in collect(Iterators.product(1:num_scales, 1:num_masses, 1:num_deltas, 1:num_Yukawas)) 
#=   
Threads.@threads for super_index in 1:num_of_parameter_points  
    inds = rolled_index(super_index, [num_scales, num_masses, num_deltas, num_Yukawas])
    i = inds[1]
    j = inds[2]
    k = inds[3]
    l = inds[4]
    =#
    Lambda_dQCD = array_scales[i]
    m_N = array_masses[j]*Lambda_dQCD
    mass_delta = array_deltas[k]
    m_L = (mass_delta+1)*m_N 
    ydark = array_Yukawas[l] #dark Yukawa coupling. No running implemented. 

    Alpha_DM = running_coupling_from_pole(m_N, Lambda_dQCD, 11*Ndark/3)
    Alpha_NN = running_coupling_from_pole(2*m_N, Lambda_dQCD, (11*Ndark-2)/3) #Alpha_DM at the scale of NN annihilation. The N flavour is active. 
    Alpha_NL = running_coupling_from_pole(m_N + m_L, Lambda_dQCD, (11*Ndark-2*2-2)/3) #Alpha_DM at the scale of NL annihilation. Both flavours are active
    Alpha_LL = running_coupling_from_pole(2*m_L, Lambda_dQCD, (11*Ndark - 2*2-2)/3) #Alpha_DM at the scale of LL annihilation
    aw_NL = running_coupling_from_scale(m_N + m_L, Mtop, alpha_W_Mtop, (19 - 4*Ndark)/6) #beta0 = 19/6 for the weak interaction below m_L
    aw_LL = running_coupling_from_scale(2*m_L, Mtop, alpha_W_Mtop, (19 - 4*Ndark)/6) #beta0 = 19/6 for the weak interaction below m_L
    aY_NL = running_coupling_from_scale(m_N + m_L, Mtop, alpha_Y_Mtop, -(41 + Ndark*4*2)/6) #beta0 = -41/6 for the U(1) interaction below m_L
    aY_LL = running_coupling_from_scale(2*m_L, Mtop, alpha_Y_Mtop, -(41 + Ndark*4*2)/6) #beta0 = -41/6 for the U(1) interaction below m_L
    BigConstant = bc_constant(m_N)

    aux1NL = 3*aw_NL + aY_NL # Auxiliary quantities which are defined to shorten the expressions for the cross sections a bit
    aux1LL = 3*aw_LL + aY_LL
    aux2 = 101*aw_LL*aw_LL + 8*aw_LL*aY_LL + 43*aY_LL*aY_LL
    aux3 = NdarkAdjoint*(NdarkAdjoint-1)

    sigmaNN = aux3/(8*Ndark*(NdarkAdjoint + 1))*pert_acs(Alpha_NN, m_N) #Derivation of these expression in the Mathematica notebook
    sigmaNL = (4*NdarkAdjoint*Alpha_NL*m_N*m_L + 2*NdarkAdjoint*Alpha_NL*m_L*m_L + m_N*m_N*(2*NdarkAdjoint*Alpha_NL + Ndark*aux1NL))/(32*pi)*pert_acs(ydark, Ndark*m_N*m_L)
    sigmaLL = 0.5*(2*pi*m_N*m_N*m_L*m_L*(Ndark*Ndark*(3*pi*aux2 - ydark*ydark*aux1LL) + 16*pi*Ndark*NdarkAdjoint*Alpha_LL*aux1LL + 8*pi*aux3*Alpha_LL*Alpha_LL) + m_L^4*(-2*pi*(NdarkAdjoint + 1)*ydark*ydark*aux1LL + 4*(NdarkAdjoint + 1)*ydark^4 + pi*pi*(16*NdarkAdjoint*Ndark*Alpha_LL*aux1LL + 8*aux3*Alpha_LL*Alpha_LL + 3*aux2*(NdarkAdjoint + 1))) + pi*pi*m_N^4*(16*Ndark*NdarkAdjoint*Alpha_LL*aux1LL + 8*aux3*Alpha_LL*Alpha_LL + 3*(NdarkAdjoint + 1)*aux2))*pert_acs(1, 8*pi*m_L*(m_N*m_N + m_L*m_L))/(Ndark*(NdarkAdjoint + 1))

    Tcrit = 0.63*Lambda_dQCD #Temperature of the phase transition
    x_PT = m_N/Tcrit 

    #Constants of FOPT (w/o factors of 1/Lambda for numerical convenience, cancel in the squeezeout step!)
    R0 = 1E-6*(Lambda_dQCD/Mpl)^(-0.9)
    R1 = (Mpl/(1E4*Lambda_dQCD))^(2/3)
    R_pocket = max(R0, R1)

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

    EquilibriumYieldN = zeros(Npoints) #Define equilibrium yields
    EquilibriumYieldL = zeros(Npoints)
    EquilibriumYieldDM = zeros(Npoints)
    for i = 1:Npoints
        EquilibriumYieldN[i] = Yeq(g_N, h_eff_dof(m_L/m_N*xvec[i]), m_L/m_N*xvec[i])
        EquilibriumYieldL[i] = Yeq(g_L, h_eff_dof(m_L/m_N*xvec[i]), m_L/m_N*xvec[i])
        EquilibriumYieldDM[i] = EquilibriumYieldN[i] + EquilibriumYieldL[i]
    end

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
    Yeq_ratios = zeros(Npoints)
    sigma_v_averaged = zeros(Npoints)
    for i in 1:Npoints
        if EquilibriumYieldN[i] > 0
            Yeq_ratios[i] = EquilibriumYieldL[i]/EquilibriumYieldN[i]
        end
        r_GS_N = 1/(1 + Yeq_ratios[i]) #ratios as defined in Griest and Seckel
        r_GS_L = Yeq_ratios[i]/(1 + Yeq_ratios[i])
        sigma_v_averaged[i] = (sigmaNN*r_GS_N*r_GS_N + 2*sigmaNL*r_GS_N*r_GS_N + sigmaLL*r_GS_N*r_GS_N) #Effective xs according to Griest and Seckel 
    end

    Wx = zeros(Npoints)
    Yx = zeros(Npoints)
    Yx[1] = EquilibriumYieldDM[1]
    Wx[1] = log(Yx[1])

    #Solution to the Boltzmann equation for the first freeze-out
    for i = 2:Npoints
        W_old = Wx[i-1]
        Wx[i] = Newton_Raphson_step(tvec[i], W_old, 0.5*BigConstant*Delta_t*sigma_v_averaged[i], g_N, h_eff_dof(m_N/xvec[i]))
    end
    Yx = exp.(Wx)

    ### FOPT: Squeezeout step ###
    Yx_squeezeout = 1.5/pi*sqrt(15*Yx[Npoints]/(2*pi*h_eff_dof(Tcrit)*R_pocket^3))

    ### Entropy dilution due to glueball decay ###
    m_glueball = 7*Lambda_dQCD #Mass of the lightest 0++ glueball

    x_freeze_out = GB_freeze_out_estimate(1, m_glueball, R_max) #freeze-out of dark gluons
    Y_GB = R_max/x_freeze_out #Relic yield of dark gluons
    T_MR = 4/3*m_glueball*Y_GB # Matter-radiation equality temperature, after which GB dominate the energy content
    x_MR = m_N/T_MR
    Alpha_DM_GB_decay = running_coupling_from_pole(m_glueball, Lambda_dQCD, 11*Ndark/3) #Dark gauge coupling at the mass scale of the GBs
    Alpha_weak_GB_decay = running_coupling_from_scale(m_glueball, Mtop, alpha_W_Mtop, 19/6) #Weak gauge coupling at the mass scale of the GBs
    Alpha_Y_GB_decay = running_coupling_from_scale(m_glueball, Mtop, alpha_Y_Mtop, -41/6) #Hypercharge gauge coupling at the mass scale of the GBs
    decay_const_GB = 3.06*m_glueball^3/(4*pi*Alpha_DM_GB_decay) #decay constant of the gluon after Juknevich. 
    Gamma_GB_dim8 = Alpha_DM_GB_decay*Alpha_DM_GB_decay/(8*pi*m_L^8)*1/3600*m_glueball^3*(decay_const_GB)^2*(Alpha_Y_GB_decay*Alpha_Y_GB_decay/4 + 1.5*Alpha_weak_GB_decay*Alpha_weak_GB_decay) #dim 8 Glueball decay rate after Juknevich.
    Gamma_GB_dim6 = (ydark*ydark*Alpha_DM_GB_decay*decay_const_GB/(m_L*m_N))^2/(72*pi^3*m_glueball)  #dim 6 Glueball decay rate after Juknevich.
    Gamma_GB = Gamma_GB_dim6 + Gamma_GB_dim8

    dil_fac = (1 + 1.65*g_average(1e-5, T_MR, 10)*cbrt(T_MR^4/(Gamma_GB*Mpl))^2)^(-0.75)
    #x_dilution = x_PT*100
    Yx_dilution = dil_fac*Yx_squeezeout

    ### Calculation of the relic abundance ###
    x_today = m_N/T0
    Omega_relic = reduced_Hubble_squared*Yx_dilution*s0*m_N/rho_crit

    write(results_file, join((m_N, Lambda_dQCD, x_PT, Alpha_DM, mass_delta, m_L, ydark, R_pocket, Yx[Npoints], Yx_squeezeout, x_freeze_out, T_MR, Gamma_GB, dil_fac, Omega_relic),","),"\n")
end

close(results_file)