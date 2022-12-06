using Plots
using LaTeXStrings

function Newton_Raphson_step(t, W_old, cs, g_deg, gS) #Method to calculate the next W=log(Y) point in the implicit backward Euler method.
    # t = log(x) is the t_(n+1) time step, W_old is W_n and cs is Delta_t*lambda(t_(n+1))*geff(t_(n+1)). g_deg is the degeneracy of the annihilating particles.
    # gS are the effective entropy dofs. 
    W_try_initial = W_old*1.
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
    90/(2*pi)^3.5*g_deg/gS*x*sqrt(x)*exp(-x)
end

#Define parameters of implicit Euler backward solution method
Delta_t = 1E-5;
x_initial = 20;
x_final = 1E2;
t_initial = log(x_initial);
t_final = log(x_final);

#Define physics parameters of the model
g_quark = 4;
m_quark = 1E4;
Alpha_DM = 0.1; 
sigma0 = pi*(Alpha_DM/m_quark)^2;

#Define constant physics parameters
Mpl = 1.221E19;
H0 = 1.447E-42;
T0 = 2.35E-13;
rho_crit = 3.724E-47;
s0 = 2.225E-38;
BigConstant = sqrt(pi/45)*Mpl*m_quark;
lambda = sigma0*BigConstant;

#First define initial conditions:
tvec = collect(t_initial:Delta_t:t_final);
xvec = exp.(tvec);
Npoints = length(tvec);

EquilibriumYield = zeros(Npoints);
for i = 1:Npoints
    EquilibriumYield[i] = Yeq(g_quark, 1, xvec[i])
end
Yx = zeros(Npoints);
Wx = zeros(Npoints);
Yx[1] = EquilibriumYield[1];
Wx[1] = log(Yx[1]);

#Solution to the Boltzmann equation
for i = 2:Npoints
    W_old = Wx[i-1];
    Wx[i] = Newton_Raphson_step(tvec[i], W_old, lambda*Delta_t, g_quark, 1);
end
Yx = exp.(Wx)

plot(xvec, [EquilibriumYield, Yx], title="WIMP freeze-out", label=[L"Y_{eq}(x)" L"Y(x)"], xlabel="x = m/T", ylabel="Y(x)", xaxis=:log, yaxis=:log, xlims = (x_initial, x_final), ylims = (1E-20, 1E-3))
savefig("FreezeOut.png")

plot(tvec, Wx)
savefig("W_plot.png")

plot(xvec, Yx, xaxis=:log, yaxis=:log)
savefig("Y_plot.png")