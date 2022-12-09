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
    #90/(2*pi)^3.5*g_deg/gS*x*sqrt(x)*exp(-x) #This is the non-relativistic version of the line below
    11.25*g_deg*(x/(pi*pi))^2/gS*besselk(2,x)
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

function sigma_v_interpolation(x)
    return sigma0
end

function aux_func(x) #Auxiliary function for the determination of xf. It gives H(xf) - Gamma(xf) ( == 0 at freeze-out)
    temp = m_quark/x
    f = BigConstant*eff_dof_sqrt(temp)*sigma_v_interpolation(x)*Yeq(g_quark, h_eff_dof(temp), x) - x
end

function freeze_out_estimate(xft) #Numerical estimate of the freeze-out x with the secant method. xft is the initial guess
    xf0 = xft #first try
    xf1 = xf0 - aux_func(xf0)*2*0.001/(aux_func(xf0+0.001)-aux_func(xf0-0.001))
    diff = abs(xf1 - xf0)
    if diff < 1E-4
        xf2 = xf1
    else
        while diff > 1E-4
            xf2 = xf1 - aux_func(xf1)*(xf1 - xf0)/(aux_func(xf1) - aux_func(xf0))
            diff = abs(xf2 - xf1)
            xf0 = copy(xf1)
            xf1 = copy(xf2)
        end
    end
    return xf2
end

#Define parameters of implicit Euler backward solution method
Delta_t = 1E-5
x_initial = 1
x_final = 1E9
t_initial = log(x_initial)
t_final = log(x_final)

#Define physics parameters of the model
g_quark = 4
m_quark = 1E4
Alpha_DM = 0.1
sigma0 = pi*(Alpha_DM/m_quark)^2

#Define constant physics parameters
Mpl = 1.221E19
H0 = 1.447E-42
T0 = 2.35E-13
rho_crit = 3.724E-47
s0 = 2.225E-38
BigConstant = sqrt(pi/45)*Mpl*m_quark

#First define initial conditions:
tvec = collect(t_initial:Delta_t:t_final)
xvec = exp.(tvec)
Npoints = length(tvec)

#Read and define degrees of freedom 
dof_file = DataFrame(CSV.File("DegreesOfFreedom.txt"))
temperatures = dof_file[!,1]
geff = dof_file[!,2]
heff = dof_file[!,3]
g_star_eff_sqrt = dof_file[!,4]
num_of_g_points = length(temperatures)
g_star_eff_vec = eff_dof_sqrt.(m_quark./xvec)

EquilibriumYield = zeros(Npoints)
for i = 1:Npoints
    EquilibriumYield[i] = Yeq(g_quark, h_eff_dof(m_quark/xvec[i]), xvec[i])
end

#Here the thermally averaged cross section is read in. 
sigma_v_averaged = ones(Npoints)
for i in 1:Npoints
    sigma_v_averaged[i] = sigma_v_interpolation(xvec[i])
end

#Estimate the freeze-out xf:
xf = freeze_out_estimate(25)
Y_infty = Yeq(g_quark, h_eff_dof(m_quark/xf), xf)
freeze_out_index = findfirst(xvec -> xvec > xf, xvec)
Y_infty_vec = zeros(Npoints)
for i = freeze_out_index:Npoints
    Y_infty_vec[i] = Y_infty
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

ytics = 10.0.^collect(range(-20, -3, step=1))
xtics = 10.0.^collect(range(log(x_initial), log(x_final), step=1))

plot(xvec, [EquilibriumYield, Yx, Y_infty_vec], title="WIMP freeze-out", label=[L"Y_{eq}(x)" L"Y(x)" L"Y_\infty = Y_{eq}(x_f)"], yticks = ytics, xticks = xtics, minorticks = 10, minorgrid = true, xlabel="x = m/T", ylabel="Y(x)", xaxis=:log, yaxis=:log, xlims = (x_initial, x_final), ylims = (1E-20, 1E-3))
savefig("FreezeOut.png")

plot(xvec, Yx, title="Relic Yield", minorticks = 10, minorgrid = true, xlabel="x = m/T", ylabel="Y(x)", yticks = ytics, xaxis=:log, yaxis=:log, xticks = xtics, xlims = (x_initial, x_final), ylims = (1E-20, 1E-3))
savefig("Y_plot.png")

