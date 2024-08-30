%% Initialization 
clc
clear all
E = zeros(3, 3); 
cc_profile = readmatrix('cc_values.xlsx');
des_profile = readmatrix('descent_values.xlsx');
%% Energy Consumption at Climb
clc 
%Input Parameters
phi = 0.19; %
beta = 20; %
V = 36; %
xh = 0;
rho = 1.056; %
R = 0.82;
d = 17000; %
omega = 200; %
r = 0.7*R;
%Calculations
phiD = rad2deg(phi);
c = 0.12*R; 
U0 = sqrt(V^2 + (omega*r)^2); 
l = V/(omega*R); 
SR = 3*c/(2*pi*r); 
AOA = round(beta - phiD); 
i = find(cc_profile(:, 1)==AOA);
k = 4; 
    CL = cc_profile(i, k+1); 
    CD = cc_profile(i, k+7); 
    f1 = 2*((CL*cos(phi) - CD*sin(phi)) * SR * (U0/V)^2); 
    f2 = @(x) 2*((CL*sin(phi) + CD*cos(phi)) * SR * (U0/V)^2 * x); 
    K_t = f1*(1-xh); 
    K_q = integral(f2, xh, 1);
    K_p = K_q/l; 
    T = K_t * (0.5*rho*pi*R^2*V^2); 
    n = K_t/K_p; 
    Ecl = (T * d)/n; 
    E(1, 2) = Ecl; 
    E(1, 3) = Ecl; 
%% Energy Consumption at Cruise
clc 
%Input Parameters
phi = 0.21; %
beta = 22; %
V = 80; %
xh = 0;
rho = 0.905; %
R = 0.82;
d = 140000; %
omega = 200; %
r = 0.7*R;
%Calculations
phiD = rad2deg(phi);
c = 0.12*R; 
U0 = sqrt(V^2 + (omega*r)^2); 
l = V/(omega*R); 
SR = 3*c/(2*pi*r); 
AOA = round(beta - phiD); 
i = find(cc_profile(:, 1)==AOA);
k = 4; 
    CL = cc_profile(i, k+1); 
    CD = cc_profile(i, k+7); 
    f1 = 2*((CL*cos(phi) - CD*sin(phi)) * SR * (U0/V)^2); 
    f2 = @(x) 2*((CL*sin(phi) + CD*cos(phi)) * SR * (U0/V)^2 * x); 
    K_t = f1*(1-xh); 
    K_q = integral(f2, xh, 1);
    K_p = K_q/l; 
    T = K_t * (0.5*rho*pi*R^2*V^2); 
    n = K_t/K_p; 
    Ecr = (T * d)/n; 
    E(2, 2) = Ecr; 
    E(2, 3) = Ecr; 
%% Energy Consumption at Descent
clc 
%Input Parameters
phi = 0.24; %
beta = 12; %
V = 30; %
xh = 0;
rho = 1.056; %
R = 0.82;
d = 26000; %
omega = 100; %
r = 0.7*R;
%Calculations
phiD = rad2deg(phi);
c = 0.12*R; 
U0 = sqrt(V^2 + (omega*r)^2); 
l = V/(omega*R); 
SR = 3*c/(2*pi*r); 
AOA = -5; 
i = find(des_profile(:, 1)==AOA);
k = 4; 
    CL = des_profile(i, k+1); 
    CD = des_profile(i, k+7); 
    f2 = @(x) 2*((CL*sin(phi) + CD*cos(phi)) * SR * (U0/V)^2 * x);  
    K_q = integral(f2, xh, 1);
    K_p = K_q/l; 
    Pw = 0.5*rho*pi*R^2*V^3;
    tw = d/V; 
    n = K_p; 
    Ew = n*Pw*tw; 
    E(3, 2) = Ew; 
k = 2; 
    CL = des_profile(i, k+1); 
    CD = des_profile(i, k+7); 
    f2 = @(x) 2*((CL*sin(phi) + CD*cos(phi)) * SR * (U0/V)^2 * x);  
    K_q = integral(f2, xh, 1);
    K_p = K_q/l; 
    Pw = 0.5*rho*pi*R^2*V^3;
    tw = d/V; 
    n = K_p; 
    Ew1 = n*Pw*tw; 
    E(3, 3) = Ew1; 
%%
E(3, 1) = 0; 
clc
names = {'No Renegeneration', 'Conventional', 'Morphing'};
bar(E', 'stacked');
set(gca, 'xticklabel', names)
legend('Climb Consumption', 'Cruise Consumption', 'Descent Consumption');
%% 
clc
phi = 0.24; %
beta = 12; %
AOA = round(beta - phiD)
%P1 = E(1, 1) + E(2, 1) 
%P2 = E(1, 2) + E(2, 2) + E(3, 2)