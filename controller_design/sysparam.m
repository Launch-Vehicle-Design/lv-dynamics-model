function param = sysparam
% file names
top_down_filename = "thumper.mat";
bottom_up_filename = "vehicle_sizing.mat";

% top down mass
load(top_down_filename,"optimal");
load(bottom_up_filename,"vehicle_sizing");
mass_mat = [optimal(1:3); optimal(6:8)];

% environmental parameter
param.mu = 3.986004418e14;
param.earthR = 6357e3;
param.OMEGA = [0;0;2*pi/(24*3600)];
param.skewOMEGA = [0 -param.OMEGA(3) 0; param.OMEGA(3) 0 0; 0 0 0];
param.gamma = 1.4;
param.g0 = 9.80665;

% launch vehicle parameter
param.S = pi*(24*0.0254)^2/4;
param.Isp1 = 278;
param.Isp2 = 369.5;
param.Isp_acs = 220;
param.OFratio = 1.275;
param.mass_mat = mass_mat;
param.vehicle_sizing = vehicle_sizing;
param = extract_bottomup(vehicle_sizing,param);

% basic thrust parameter
param.TtoW_1st = 1.6;
param.TtoW_2nd = 0.8;
param.T_acs = 5;
param.maxT_1st = param.TtoW_1st*param.m0*param.g0;
param.maxT_2nd = param.TtoW_2nd*param.m02*param.g0;
param.T_limit = [100 100 60 100]/100;

% TVC pointing matrix
param.tvc_yaw = @(yaw) [cos(yaw) sin(yaw) 0; -sin(yaw) cos(yaw) 0; 0 0 1];
param.tvc_pitch = @(pitch) [cos(pitch) 0 -sin(pitch); 0 1 0; sin(pitch) 0 cos(pitch)];
param.tvc_point = @(y,p) [cos(p).*cos(y); cos(p).*sin(y); -sin(p)];
param.tvc_limit = 7/180*pi*[-1 1];
% LITVC flow rate relations
relation_table = [
    0 0.02 0.03 0.05 0.06 0.075 0.08 0.09;
    0 0.05 0.075 0.15 0.21 0.3 0.33 0.4];
param.litvc_poly = polyfit(relation_table(1,:),relation_table(2,:),2);
param.litvc_mdot = @(fsfa) polyval(param.litvc_poly,fsfa);
param.litvc_limit = 3/180*pi*[-1 1];

% ACS thruster pointing
param.acs_yz_point = [11 9 15 13 19 17 23 21]/8*pi;
param.acs_xy_point = pi/4;
param.acs_point = [cos(param.acs_xy_point)*ones(size(param.acs_yz_point));
    sin(param.acs_xy_point)*cos(param.acs_yz_point);
    sin(param.acs_xy_point)*sin(param.acs_yz_point)];
param.r_acs_mount = [-param.d.enge2nd*ones([1,8]); cos(pi/8:pi/4:2*pi); sin(pi/8:pi/4:2*pi)];

% CP estimation by Barrowman's equation
param = cp_est(param);

% Grid fin parameters
param.dclda = @(M) min(4./sqrt(M.^2-1),2*pi);
param.gfS = param.brm.S*param.brm.C_R;
param.gfn = param.brm.N;
param.gf_point = [zeros(1,4); cos([3*pi/4 pi/4 3*pi/4 pi/4]); sin([3*pi/4 pi/4 3*pi/4 pi/4])];
param.r_gf_mount = [-param.d.enge1st*ones([1,4]); cos(pi/4:pi/2:2*pi); sin(pi/4:pi/2:2*pi)];
param.gf_limit = [-45 45]*pi/180;

% initial release condition
param.release.lat = 80/180*pi;
param.release.long = -180/180*pi;
param.release.alti = 40000*0.3048;
param.release.velo = 250.786;

% final condition parameter
param.final_ref.v = 7754.1;
param.final_ref.alt = 251460;
param.final_ref.fpa = 0;
param.final_ref.orb_inc = 100/180*pi;
param.final_ref.plf_dropalti = 140e3;

% body shape
[upplf,botplf] = xcone(12*0.0254,param.l.plf,[0;0;0],1);
body_len = param.vehicle_sizing.lengths(end)-param.l.plf;
[up1st,bot1st,side1st] = xtube(12*0.0254,body_len,[-param.vehicle_sizing.lengths(end);0;0]);
param.body.upplf = upplf;
param.body.botplf = botplf;
param.body.up1st = up1st;
param.body.bot1st = bot1st;
param.body.side1st = side1st;

% vehicle graph plotter test
% ang = pi/4;
% C = [cos(ang) 0 -sin(ang); 0 1 0; sin(ang) 0 cos(ang)];
% lvgraph(lvrot(C,param.body));
end

