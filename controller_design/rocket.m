%% FUNCTION - Dynamics of a launch vehicle
function [dxdt,forces] = rocket(t,x,u,param,phase)
% % x represents states as follows 7x1
% 1:3 - r, position of the launch vehicle, iner
% 4:6 - v, velocity of the launch vehicle, iner
% 7 - m, mass of the launch vehicle

% % u represents control as follows 4x1
% 1:3 - TVC direction, iner
% 4 - throttle

r = x(1:3,:); v = x(4:6,:); m = x(7,:);
tvc = u(1:3,:); throtl = u(4,:);

% norm array
rn = vecnorm(r);

% derived parameters
h = rn-param.earthR;
atmo_profile = atmo(h*param.scales.length);
rho = atmo_profile.rho/param.scales.density;
a = sqrt(param.gamma.*atmo_profile.P./rho/param.scales.pressure);

% aerodynamic drag
atmo_v = param.skewOMEGA*r;
vrel = v - atmo_v; vreln = vecnorm(vrel);
e_d = -vrel./vreln;
mach = vreln./a;
Cd = 0.5*ones(size(mach));
% for i = 1:length(mach) Cd(i) = CD(mach(i),0); end

% external forces
ag = -param.mu./rn.^3.*r;
Fd = 1/2.*rho.*vreln.^2.*Cd*param.S.*e_d;
Ft = zeros(size(ag)); mdot = zeros(size(rn));
if phase == 1
    Ft = throtl.*param.maxT_1st.*tvc;
    mdot = vecnorm(Ft)./param.Isp1./param.g0;
elseif phase == 2
    Ft = throtl.*param.maxT_2nd.*tvc;
    mdot = vecnorm(Ft)./param.Isp2./param.g0;
end
% total force
atot = ag + (Fd + Ft)./m;

% differential of states
dxdt = zeros(size(x));
dxdt(1:3,:) = v;
dxdt(4:6,:) = atot;
dxdt(7,:) = -mdot;

forces = [ag; Fd; Ft];
