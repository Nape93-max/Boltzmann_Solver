function fRK4(yn, xn)
    K1 = -g*xn;
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
Delta_x = 1E-1;
x_initial = 0.;
x_final = 4.;
Npoints = ceil(Int8,(x_final-x_initial)/Delta_x);

#First define initial conditions:
g = 9.81;
Y0 = 100*ones(Npoints);
Yx = copy(Y0);
xvec = collect(x_initial:Delta_x:x_final);

for i in 2:Npoints
    Yx[i] = RK4step(Yx[i-1], xvec[i-1]);
end
