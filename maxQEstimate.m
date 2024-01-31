clear; clc; close all
addpath("parameter_functions\")

start_burn_alt = 39599*0.3048;
acceleration = 1.6*9.80665;
start_vel = 220;
mass_init = 2250;
isp = 270;

avg_climb_angle = 35;
aca_rad = avg_climb_angle/180*pi;

gamma = 1.4;

t = 0:100;
m(1) = mass_init; v(1) = start_vel; h(1) = start_burn_alt;
dt = diff(t);
for ind = 2:length(t)
    acce = acceleration - g(h(ind-1))*sin(aca_rad) - 0.2*v(ind-1)^2*0.5*rho(h(ind-1))*(pi*0.3048^2)/m(ind-1);
    v(ind) = v(ind-1) + acce*dt(ind-1);
    h(ind) = h(ind-1) + (v(ind-1)*dt(ind-1) + 1/2*acce*dt(ind-1)^2)*sin(aca_rad);
    m(ind) = m(ind-1) - acceleration*m(ind-1)/(isp*g(h(ind-1)));
end
% v = start_vel + acceleration*t;
% h = start_burn_alt + (start_vel*t + 1/2*acceleration*t.^2)*sin(aca_rad);

ps = P(h);
a = sqrt(gamma*ps./rho(h));

Q = ps.*(1+(gamma-1)/2*(v./a).^2).^(gamma/(gamma-1))-ps;
[val, ind] = max(Q);

disp("Max Q estimated value (psi): "+num2str(max(Q)/6895))
disp("Max Q estimated time (s): "+num2str(t(ind)))
disp("Max Q estimated velocity (ft/s): "+num2str(v(ind)/0.3048))
disp("Max Q estimated altitude (ft): "+num2str(h(ind)/0.3048))

disp("Max Q estimated velocity (m/s): "+num2str(v(ind)))
disp("Max Q estimated altitude (m): "+num2str(h(ind)))
