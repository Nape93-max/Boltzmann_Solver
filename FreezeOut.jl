using SpecialFunctions
using CSV, DataFrames

#Define constant physics parameters
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

#BEGIN function definitions:

function Newton_Raphson_step(x, W_old, cs, g_deg, gS) #Method to calculate the next W=log(Y) point in the implicit backward Euler method.
    # x is the x_(n+1) time step, W_old is W_n and cs is Delta_t*lambda(t_(n+1))*geff(t_(n+1)). g_deg is the degeneracy of the annihilating particles.
    # gS are the effective entropy dofs. 
    W_try_initial = W_old
    W_new = Newton_Raphson_iteration(x, W_old, cs, W_try_initial, g_deg, gS)
    diff = abs(log(abs(W_new/W_try_initial)))
    while diff > 1E-2 
        W_try = deepcopy(W_new);
        W_new = Newton_Raphson_iteration(x, W_old, cs, W_try, g_deg, gS)
        diff = abs(log(abs(W_new/W_try)))
    end
    if isnan(W_new)
        println("ALARM: W_new in function Newton_Raphson is NaN. Maybe relax bound on diff in the while loop as a fix.")
    end
    return W_new
end

function Newton_Raphson_iteration(x, W_old, cs, W_previous, g_deg, gS) #Does one NR-step to calculate a trial W^(i+1)_(n+1). W_previous = W^(i)_(n+1)
    A = exp(W_previous);
    B = Yeq(g_deg, gS, x)^2/A 
    W_next = W_previous - (W_previous - W_old + cs*(A - B))/(1 + cs*(A + B))
    return W_next
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

function aux_func(x, m, acs, cs_xvalues, cs_y_values, coeffs) #Auxiliary function for the determination of xf. It gives H(xf) - Gamma(xf) ( == 0 at freeze-out)
    temp = m/x
    f = BigConstant*eff_dof_sqrt(temp)*sigma_v_cspline(x, cs_xvalues, cs_y_values, acs, coeffs)*Yeq(g_quark, h_eff_dof(temp), x) - x
end

function freeze_out_estimate(xft, m, acs, cs_xvalues, cs_yvalues, coeffs) #Numerical estimate of the freeze-out x with the secant method. xft is the initial guess, m is the DM mass and acs is the constant annihilation xs. cs_xvalues are the x_values of the fully dressed acs. cs_yvalues are the corresponding acs values. coeffs are the spline interpolation coefficients.
    xf0 = xft #first try
    xf1 = xf0 - aux_func(xf0, m, acs, cs_xvalues, cs_yvalues, coeffs)*2*0.001/(aux_func(xf0+0.001, m, acs, cs_xvalues, cs_yvalues, coeffs)-aux_func(xf0-0.001, m, acs, cs_xvalues, cs_yvalues, coeffs))
    diff = abs(xf1 - xf0)
    if diff < 1E-4
        xf2 = xf1
    else
        while diff > 1E-4
            xf2 = xf1 - aux_func(xf1, m, acs, cs_xvalues, cs_yvalues, coeffs)*(xf1 - xf0)/(aux_func(xf1, m, acs, cs_xvalues, cs_yvalues, coeffs) - aux_func(xf0, m, acs, cs_xvalues, cs_yvalues, coeffs))
            diff = abs(xf2 - xf1)
            xf0 = copy(xf1)
            xf1 = copy(xf2)
        end
    end
    return xf2
end

function pert_acs(alpha, m)
    h = alpha/m
    return pi*h*h
end

function bc_constant(m)
    return sqrt(pi/45)*Mpl*m
end

function sigma_v_interpolation(x_input, y_input) # Interpolation of sigma_v with a cubic spline.
    # x(y)_input are the vectors containing input data of the xs. cs is a constant part of the xs.
    # The output is a 3*n array of (beta; gamma; delta) coefficients for the spline.
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

function sigma_v_cspline(x, x_input, y_input, cxs, coeffs) #Cubic spline interpolation step with the previously calculated coefficient vectors beta, gamma and delta
    # encoded in coeffs. x is the argument, x(y)_input is the input x(y)-data and cxs is the constant part of the cross section
    n = length(x_input)
    xvals = log.(x_input) #logarithmic interpolation
    yvals = log.(y_input)
    if x >= last(x_input)
        ind = n
    else
        ind = findfirst(x_input -> x_input > x, x_input) #find relevant interval
    end

    beta = coeffs[ind] #coefficients beta_i, gamma_i and delta_i
    gamma = coeffs[ind + n]
    delta = coeffs[ind + 2*n]

    arg = log(x) - xvals[ind] # x-x_i
    return cxs*exp(yvals[ind] + arg*(beta + arg*(gamma + arg*delta))) #Cubic spline step
end

function running_coupling_from_pole(Q, Lambda, beta0) #Running coupling at the scale Q given a Landau pole Lambda with a theory with a beta0 coefficient
    return 2*pi/(beta0*log(Q/Lambda))
end

function running_coupling_from_scale(Q, mu, alphamu, beta0) #Running coupling at the scale Q given the running alphamu at another scale mu with a beta0 coefficient. 
    return 1/(1/alphamu + beta0*log(Q/mu)/(2*pi))
end

function entropy_density(T) #Returns the entropy density of a given species
    s = 2*pi*pi/45*h_eff_dof(T)*T^3
end

function quark_freeze_out(x_final, m, sigma, BC, g_dm)
    #Define parameters of implicit Euler backward solution method
    Delta_t = 1E-4
    x_initial = 1E-1
    t_initial = log(x_initial)
    t_final = log(x_final)

    #First define initial conditions:
    tvec = collect(t_initial:Delta_t:t_final)
    xvec = exp.(tvec)
    Npoints = length(tvec)

    g_star_eff_vec = eff_dof_sqrt.(m./xvec) #Sqrt of effective degrees of freedom
    h_eff_dof_vec = h_eff_dof.(m./xvec) #Effective entropic dofs

    EquilibriumYield = zeros(Npoints)
    for i = 1:Npoints
        EquilibriumYield[i] = Yeq(g_dm, h_eff_dof_vec[i], xvec[i])
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
        sigma_v_averaged[i] = sigma
    end

    Wx = zeros(Npoints)
    Yx = zeros(Npoints)
    Yx[1] = EquilibriumYield[1]
    Wx[1] = log(Yx[1])

    modified_BC = BC*Delta_t*0.5 #Modified BC make code more efficient    
    #Solution to the Boltzmann equation for the first freeze-out
    for i = 2:Npoints
        W_old = Wx[i-1]
        Wx[i] = Newton_Raphson_step(xvec[i], W_old, modified_BC*g_star_eff_vec[i]/xvec[i]*sigma_v_averaged[i], g_dm, h_eff_dof_vec[i])
    end

    Yx = exp.(Wx)
    return Yx[Npoints-1]
end
