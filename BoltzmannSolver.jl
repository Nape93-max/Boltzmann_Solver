function fRK4(yn, xn)
    K1 = lambda*exp(-xn)*(Yeq(xn)^2*exp(-yn)-exp(yn));
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

function Yeq(x)
    y = 0.145*x^1.5*exp(-x)
end

#Define parameters of RK4 method
Delta_x = 1E-4;
x_initial = 2.;
x_final = 10;
Npoints = ceil(Int64,(x_final-x_initial)/Delta_x);

#Define physics parameters
lambda = 1E10;

#Define initial conditions:
Yx = zeros(Npoints);
xvec = collect(x_initial:Delta_x:x_final);
Yx[1] = log(Yeq(xvec[1]));

for i in 2:Npoints
    Yx[i] = RK4step(Yx[i-1], xvec[i-1]);
end

Yx