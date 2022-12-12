using Plots
using LaTeXStrings
using SpecialFunctions
using CSV, DataFrames

function Newton_Raphson_step(t, W_old, cs, g_deg, gS) #Method to calculate the next W=log(Y) point in the implicit backward Euler method.
    # t = log(x) is the t_(n+1) time step, W_old is W_n and cs is Delta_t*lambda(t_(n+1))*geff(t_(n+1)). g_deg is the degeneracy of the annihilating particles.
    # gS are the effective entropy dofs. 
    W_try_initial = W_old
    W_new = Newton_Raphson_iteration(t, W_old, cs, W_try_initial, g_deg, gS)
    diff = abs(log(W_new/W_try_initial))
    while diff > 1E-4 
        W_try = copy(W_new);
        W_new = Newton_Raphson_iteration(t, W_old, cs, W_try, g_deg, gS)
        diff = abs(log(abs(W_new/W_try)))
    end
    return W_new
end

function Newton_Raphson_iteration(t, W_old, cs, W_previous, g_deg, gS) #Does one NR-step to calculate a trial W^(i+1)_(n+1). W_previous = W^(i)_n
    A = exp(W_previous);
    B = Yeq(g_deg, gS, exp(t))^2/A
    W_next = W_previous - (W_old - W_previous + cs*(A - B))/(1 + cs*(A + B))
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

function aux_func(x, m, acs, cs_xvalues, cs_y_values) #Auxiliary function for the determination of xf. It gives H(xf) - Gamma(xf) ( == 0 at freeze-out)
    temp = m/x
    f = BigConstant*eff_dof_sqrt(temp)*sigma_v_interpolation(acs, x, cs_xvalues, cs_y_values)*Yeq(g_quark, h_eff_dof(temp), x) - x
end

function freeze_out_estimate(xft, m, acs, cs_xvalues, cs_yvalues) #Numerical estimate of the freeze-out x with the secant method. xft is the initial guess, m is the DM mass and acs is the constant annihilation xs. cs_xvalues are the x_values of the fully dressed acs. cs_yvalues are the corresponding acs values.
    xf0 = xft #first try
    xf1 = xf0 - aux_func(xf0, m, acs, cs_xvalues, cs_yvalues)*2*0.001/(aux_func(xf0+0.001, m, acs, cs_xvalues, cs_yvalues)-aux_func(xf0-0.001, m, acs, cs_xvalues, cs_yvalues))
    diff = abs(xf1 - xf0)
    if diff < 1E-4
        xf2 = xf1
    else
        while diff > 1E-4
            xf2 = xf1 - aux_func(xf1, m, acs, cs_xvalues, cs_yvalues)*(xf1 - xf0)/(aux_func(xf1, m, acs, cs_xvalues, cs_yvalues) - aux_func(xf0, m, acs, cs_xvalues, cs_yvalues))
            diff = abs(xf2 - xf1)
            xf0 = copy(xf1)
            xf1 = copy(xf2)
        end
    end
    return xf2
end

function pert_acs(alpha, m)
    h = alpha/m
    sigma = pi*h*h
end

function bc_constant(m)
    c = sqrt(pi/45)*Mpl*m
end

#Define parameters of implicit Euler backward solution method
Delta_t = 1E-4
x_initial = 1E-1
x_final = 1E5
t_initial = log(x_initial)
t_final = log(x_final)

#Define physics parameters of the model
g_quark = 4
m_quark = 1E4
Alpha_DM = 0.1
sigma0 = pert_acs(Alpha_DM, m_quark)

#Define constant physics parameters
const Mpl = 1.221E19
const H0 = 1.447E-42
const T0 = 2.35E-13
const rho_crit = 3.724E-47
const s0 = 2.225E-38
BigConstant = bc_constant(m_quark)

#First define initial conditions:
tvec = collect(t_initial:Delta_t:t_final)
xvec = exp.(tvec)
Npoints = length(tvec)

#Read and define degrees of freedom 
dof_file = DataFrame(CSV.File("DegreesOfFreedom.txt"))
const temperatures = dof_file[!,1]
const geff = dof_file[!,2]
const heff = dof_file[!,3]
const g_star_eff_sqrt = dof_file[!,4]
const num_of_g_points = length(temperatures)
g_star_eff_vec = eff_dof_sqrt.(m_quark./xvec)

EquilibriumYield = zeros(Npoints)
for i = 1:Npoints
    EquilibriumYield[i] = Yeq(g_quark, h_eff_dof(m_quark/xvec[i]), xvec[i])
end

function sigma_v_interpolation(cs, x_input, y_input) # Interpolation of sigma_v with a cubic spline. x(y)_input are the vectors containing input data of the xs. cs is a constant part of the xs. The output is a 3*n array of (beta; gamma; delta) coefficients for the spline.
    n = length(x_input) #read in data in logarithmic scaling for improved accuracy
    xvals = log.(x_input)
    yvals = log.(y_input)

    beta = zeros(n) # initialise output vectors
    gamma = zeros(n)
    delta = zeros(n)

    A = zeros(n) # Initialise auxiliary vectors for the spline
    M = zeros(n)

    h = prepend!(diff(xvals), xvals[1] - xvals[n]) #Set matrix entries for the spline matrix problem
    hfirst = h[1] #auxiliary variable
    y_diff_vec = prepend!(diff(yvals), yvals[1] - yvals[n]) #Auxiliary expr appearing in the definition of teh inhomogeneity d_i
    lambda = ones(n) 
    d = zeros(n)
    mu = zeros(n)

    for i = 1:n-1
        lambda[i] = 1/(1 + h[i]/h[i+1])
        mu[i] = 1 - lambda[i] 
        d[i] = 6/(h[i] + h[i+1])*((yvals[i+1] - yvals[i])/h[i+1] - y_diff_vec[i]/h[i])
    end
    lambda[n] = 1/(1 + h[n]/h[1])
    mu[n] = 1 - lambda[n]
    d[n] = 6/(h[1] + h[n])*((yvals[1] - yvals[n])/h[1] - y_diff_vec[n]/h[n])

    ### Solution via Gauss elimination
    
    ###
    
    M_diff_vec = append!(diff(M), M[1] - M[n])

    for i = 1:n-1 #Now calculate the coefficients beta_i, gamma_i and delta_i from the solution of the matrix problem
        hind = h[i+1]
        A[i] = y_diff_vec[i+1]/hind - hind/6*M_diff_vec[i]
        beta[i] = A[i] - M[i]*hind*hind/6
        delta[i] = M_diff_vec[i]/(6*hind)
    end
    A[n] = y_diff_vec[1]/hfirst - hfirst/6*M_diff_vec[n] #Declare the boundary terms
    beta[n] = yvals[n] - M[n]*hfirst*hfirst/6
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

#Here the thermally averaged cross section is read in. 
sigma_v_file = DataFrame(CSV.File("Sigma_eff.txt"))
sigma_v_x_values = sigma_v_file[!,1]
sigma_v_y_values = sigma_v_file[!,2]

sigma_v_averaged = ones(Npoints) #Initialise array for interpolated annihaltion xs
sigma_v_averaged_coeffs = sigma_v_interpolation(sigma0, sigma_v_x_values, sigma_v_y_values) #Cubic spline fit coefficient vector (beta, gamma, delta)

for i in 1:Npoints
    sigma_v_averaged[i] = sigma_v_cspline(xvec[i], sigma_v_x_values, sigma_v_y_values, sigma0, sigma_v_averaged_coeffs)
end

plot([xvec, sigma_v_x_values], [sigma_v_averaged, sigma0.*sigma_v_y_values], xaxis=:log, yaxis=:log)
savefig("sigma_v_data.png")

#=
#Estimate the freeze-out xf:
xf = freeze_out_estimate(25, m_quark, sigma0, sigma_v_x_values, sigma_v_y_values)
Y_infty = Yeq(g_quark, h_eff_dof(m_quark/xf), xf)
#Y_infty2 = 3.79*xf/(eff_dof_sqrt(m_quark/xf)*Mpl*m_quark*sigma0) # This is eq. (5.45) in Kolb and Turner. Coincides with Y_infty. 
freeze_out_index = findfirst(xvec -> xvec > xf, xvec)
Y_infty_vec = zeros(Npoints)
#Y_infty2_vec = zeros(Npoints)
for i = freeze_out_index:Npoints
    Y_infty_vec[i] = Y_infty
    #Y_infty2_vec[i] = Y_infty2
end

Wx = zeros(Npoints)
Yx = zeros(Npoints)
Yx[1] = EquilibriumYield[1]
Wx[1] = log(Yx[1])

#Solution to the Boltzmann equation
for i = 2:Npoints
    W_old = Wx[i-1]
    Wx[i] = Newton_Raphson_step(tvec[i], W_old, BigConstant*Delta_t*sigma_v_averaged[i], g_quark, h_eff_dof(m_quark/xvec[i]))
end
Yx = exp.(Wx)

### Here the plotting business starts 

ytics = 10.0.^collect(range(-20, -1, step=1))
xtics = 10.0.^collect(range(log10(x_initial), log10(x_final), step=1))

plot(sigma_v_x_values, sigma_v_y_values, xaxis=:log, yaxis=:log)
savefig("sigma_v_data.png")

plot(xvec, [EquilibriumYield, Yx, Y_infty_vec], title="WIMP freeze-out", label=[L"Y_{eq}(x)" L"Y(x)" L"Y_\infty = Y_{eq}(x_f)"], yticks = ytics, xticks = xtics, minorticks = 10, minorgrid = true, xlabel="x = m/T", ylabel="Y(x)", xaxis=:log, yaxis=:log, xlims = (x_initial, x_final), ylims = (1E-16, 1E-1))
plot!([xf], seriestype = :vline, label = L"x_f")
savefig("FreezeOut.png")

plot(xvec, Yx, title="Relic Yield", minorticks = 10, minorgrid = true, xlabel=L"x = m/T", ylabel=L"Y(x)", label = L"Y(x)", yticks = ytics, xaxis=:log, yaxis=:log, xticks = xtics, xlims = (x_initial, x_final), ylims = (1E-16, 1E-1))
plot!([xf], seriestype = :vline, label = L"x_f")
savefig("Y_plot.png")

=#