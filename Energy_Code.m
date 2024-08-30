%% Initialization 
clc
clear all
E = zeros(3, 6); 
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
r = R;
%Calculations
phiD = rad2deg(phi);
c = 0.05*R; 
U0 = sqrt(V^2 + (omega*r)^2); 
l = V/(omega*R); 
SR = 3*c/(2*pi*r); 
AOA = round(beta - phiD); 
i = find(cc_profile(:, 1)==AOA);
k = 1; 
while k<=6
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
    E(1, k) = Ecl; 
    k = k+1;
end 
%% Energy Consumption at Cruise
clc 
%Input Parameters
phi = 0.21; %
beta = 20; %
V = 80; %
xh = 0;
rho = 0.905; %
R = 0.82;
d = 140000; %
omega = 200; %
r = R;
%Calculations
phiD = rad2deg(phi);
c = 0.05*R; 
U0 = sqrt(V^2 + (omega*r)^2); 
l = V/(omega*R); 
SR = 3*c/(2*pi*r); 
AOA = round(beta - phiD); 
i = find(cc_profile(:, 1)==AOA);
k = 1; 
while k<=6
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
    E(2, k) = Ecl; 
    k = k+1;
end 
%% Energy Consumption at Descent
clc 
%Input Parameters
phi = 0.25; %
beta = 20; %
V = 30; %
xh = 0;
rho = 1.056; %
R = 0.82;
d = 26000; %
omega = 100; %
r = R;
%Calculations
c = 0.05*R; 
U0 = sqrt(V^2 + (omega*r)^2); 
l = V/(omega*R); 
SR = 3*c/(2*pi*r); 
AOA = round(beta - phiD); 
i = find(des_profile(:, 1)==AOA);
k = 1; 
while k<=6
    CL = des_profile(i, k+1); 
    CD = des_profile(i, k+7); 
    f2 = @(x) 2*((CL*sin(phi) + CD*cos(phi)) * SR * (U0/V)^2 * x);  
    K_q = integral(f2, xh, 1);
    K_p = K_q/l; 
    Pw = 0.5*rho*pi*R^2*V^3;
    tw = d/V; 
    n = K_p; 
    Ew = n*Pw*tw; 
    E(3, k) = Ew; 
    k = k+1;
end 
%%
clc
names = {'0012', '0015', '4412', '4415', '2412', '2415'};
bar(E', 'stacked');
set(gca, 'xticklabel', names)
legend('Climb Consumption', 'Cruise Consumption', 'Descent Consumption');
%% 