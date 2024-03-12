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
param.Isp1 = 293.49;
param.Isp2 = 369.5;
param.OFratio = 1.325;
param.mass_mat = mass_mat;
param.vehicle_sizing = vehicle_sizing;
param = extract_bottomup(vehicle_sizing,param);

% basic thrust parameter
param.TtoW_1st = 1.6; 
param.TtoW_2nd = 0.8;
param.maxT_1st = param.TtoW_1st*param.m0*param.g0;
param.maxT_2nd = param.TtoW_2nd*param.m02*param.g0;

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
end

