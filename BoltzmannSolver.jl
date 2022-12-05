function fRK4(yn, xn)
    K1 = -exp(-xn)/25*(1E10*Yeq*exp(-25*exp(xn))*sinh(yn)+1);
end

function RK4K2(yn, xn)
    h = Delta_x/2;
    K1 = fRK4(yn, xn);
    K2 = fRK4(yn + h*K1, xn + h);
end

function RK4K3(yn, xn)
    h = Delta_x/2;
    K2 = RK4K2(yn, xn);
    K3 = fRK4(yn + h*K2, xn + h);
end

function RK4K4(yn, xn)
    K3 = RK4K3(yn, xn);
    fRK4(yn + K3*Delta_x, xn + Delta_x);
end 

function RK4step(yn, xn)
    K1 = fRK4(yn, xn);
    K2 = RK4K2(yn, xn);
    K3 = RK4K3(yn, xn);
    K4 = RK4K4(yn, xn);
    ynew = yn + Delta_x*(K1 + 2*K2 + 2*K3 + K4)/6;
end

#Define parameters of RK4 method
Delta_x = 1E-6;
x_initial = -3.;
x_final = 3;
Npoints = ceil(Int64,(x_final-x_initial)/Delta_x);

#Define physics parameters of the model
g_quark = 4;
m_quark = 1E4;
Alpha_DM = 0.1; 
sigma0 = pi*(Alpha_DM/m_quark)^2;
Yeq = 0.1;

#Define constant physics parameters
Mpl = 1.221E19;
H0 = 1.447E-42;
T0 = 2.35E-13;
rho_crit = 3.724E-47;
s0 = 2.225E-38;
BigConstant = sqrt(pi/45)*Mpl*m_quark;

#First define initial conditions:
Y0 = Yeq*ones(Npoints);
Yx = zeros(Npoints);
xvec = collect(x_initial:Delta_x:x_final);

for i in 2:Npoints
    Yx[i] = RK4step(Yx[i-1], xvec[i-1]);
end

Yx