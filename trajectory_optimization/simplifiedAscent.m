function [dxdt,forces] = simplifiedAscent(t, x, stage, firing, param, envir)
% x represents stats as follows
% 1 - v, velocity
% 2 - gamma, flight path angle
% 3 - h, altitude
% 4 - x, downrange
% 5 - m, total mass

if ~exist("stage","var")
    stage = 1;
    Isp = param.Isp1;
    m0 = param.m0;
    tbyw = param.tbyw1;
elseif stage == 1
    Isp = param.Isp1;
    m0 = param.m0;
    tbyw = param.tbyw1;
elseif stage == 2
    Isp = param.Isp2;
    m0 = param.m02;
    tbyw = param.tbyw2;
end

h = x(3); v = x(1); Rp = envir.Rp;
density = rho(h); q = 1/2*density*v^2;
p = P(h); gh = g(h); mach = v/sqrt(1.4*p/density); 
Cd = CD(mach,0); D = q*param.S*Cd;
dxdt = zeros(size(x));

if ~firing
    T = 0;
else
    T = m0*envir.g0*tbyw;
end

dxdt(1) = (T-D)/x(5)-gh*sin(x(2));
dxdt(2) = -(gh/v-v/(Rp+h))*cos(x(2));
dxdt(3) = x(1)*sin(x(2));
dxdt(4) = x(1)*cos(x(2))/(1+h/Rp);
dxdt(5) = -T/envir.g0/Isp;

forces = [T D];