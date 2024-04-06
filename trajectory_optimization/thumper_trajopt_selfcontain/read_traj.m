clear; clc; close all
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

file_name = "thumper_straj_cg3dof.mat";
if exist(file_name,"file")
    load(file_name,"log_x","log_param");
end

% set mission time since the drop
t = 0:1:600;
states = nan([log_param.nstate,length(t)]);
ctrls = nan([log_param.nctrl,length(t)]);
forces = nan([9,length(t)]);
stages = nan([1,length(t)]);

for i = 1:length(t)
    cont_traj = cont3d_traj_interp(t(i),log_x,log_param);
    states(:,i) = cont_traj.state; ctrls(:,i) = cont_traj.ctrl;
    forces(:,i) = cont_traj.force; stages(i) = cont_traj.stage;
end

earth_radius = log_param.earthR*log_param.scales.length;
rn = vecnorm(states(1:3,:)); vn = vecnorm(states(4:6,:)); alti = rn-earth_radius;
tvc = ctrls(1:3,:); throtl = ctrls(4,:); pa = acos(dot(states(4:6,:)./vn,tvc));

% body frame unit vector
vec_n = cross(states(1:3,:),states(4:6,:)); e_nx = vec_n./vecnorm(vec_n);
e_vz = states(4:6,:)./vn; e_upy = cross(e_vz,e_nx);

force_acce = forces(4:6,:)./states(7,:) + forces(7:9,:)./states(7,:);
for i = 1:size(e_nx,2)
    b_mat_b2i = [e_nx(:,i) e_upy(:,i) e_vz(:,i)];
    force_acce_b(:,i) = b_mat_b2i\force_acce(:,i);
    acce_b(:,i) = b_mat_b2i\(forces(1:3,i)+force_acce(:,i));
end

acce = forces(1:3,:) + force_acce;
gacce_vec_b = acce_b/(log_param.g0*log_param.scales.acceleration);
gs_load_vec_b = force_acce_b/(log_param.g0*log_param.scales.acceleration);

gacce_lateral = sqrt(gacce_vec_b(1,:).^2+gacce_vec_b(2,:).^2);
gacce_axial = gacce_vec_b(3,:); gs_acce = vecnorm(gacce_vec_b);
gload_lateral = sqrt(gs_load_vec_b(1,:).^2+gs_load_vec_b(2,:).^2);
gload_axial = gs_load_vec_b(3,:); gs_load = vecnorm(gs_load_vec_b);

% predict percise elliptical orbit
rf = states(1:3,end); vf = states(4:6,end);
elements = orbelm(rf,vf,log_param);
rot3 = @(a) [[cos(a) -sin(a) 0];[sin(a) cos(a) 0];[0 0 1]];
rot1 = @(a) [[1 0 0];[0 cos(a) -sin(a)];[0 sin(a) cos(a)]];
orb_rotmat = rot3(elements.raan*pi/180)*rot1(elements.inc*pi/180)*rot3(elements.omega*pi/180);
nu = (0:360)/180*pi; p = elements.a*(1-elements.e^2);
rorb = @(nuorb) p./(1+elements.e.*cos(nuorb)).*[cos(nuorb);sin(nuorb);zeros(size(nuorb))];
vorb = @(nuorb) sqrt(elements.mu/p)*[-sin(nuorb);elements.e+cos(nuorb);zeros(size(nuorb))];
rorb_rot = orb_rotmat*rorb(nu); vorb_rot = orb_rotmat*vorb(nu);

% translate to LLA
lla_rot = ecef2lla(rorb_rot',0,earth_radius)';
lla_traj = ecef2lla(states(1:3,:)',0,earth_radius)';

figure; plot3(rorb_rot(1,:),rorb_rot(2,:),rorb_rot(3,:),"Color",[0.5 0.5 0.5],"LineWidth",1,"LineStyle","--");
hold on; [X,Y,Z] = sphere(100); xl = xlim; yl = ylim; zl = zlim;
plot3(states(1,:),states(2,:),states(3,:),"Color",[1 0 0],"LineWidth",1.5);
surf(X*log_param.scales.length,Y*log_param.scales.length,Z*log_param.scales.length, ...
    "EdgeColor","none","FaceColor",[0.8 0.8 0.8],"FaceAlpha",1);
axis equal; xlim(xl); ylim(yl); zlim(zl);
xlabel("x"); ylabel("y"); zlabel("z");

uif = uifigure; g = geoglobe(uif);
geoplot3(g,lla_rot(1,:),lla_rot(2,:),lla_rot(3,:),"LineWidth",2,"Color",[0 104 58]/255); hold(g,'on')
geoplot3(g,lla_traj(1,:),lla_traj(2,:),lla_traj(3,:),"LineWidth",2,"Color",[1 1 0]);

