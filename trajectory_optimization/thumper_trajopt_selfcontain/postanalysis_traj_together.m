clear; clc; close all
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',12)

file_name = "thumper_straj_cg3dof.mat";
if exist(file_name,"file")
    load(file_name,"log_x","log_param");
end

% set mission time since the drop
dt = 0.02;
t = 0:dt:510;
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
tvc = ctrls(1:3,:); throtl = ctrls(4,:); 
tvca = acos(dot(states(4:6,:),tvc)./vn./vecnorm(tvc));
pa = pi/2-acos(dot(states(4:6,:),states(1:3,:))./vn./rn);

% downrange
init_pos = states(1:3,1);
dr = earth_radius*acos(init_pos'*states(1:3,:)/norm(init_pos)./rn);

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

% dynamic pressure
atmo_profile = atmo(alti);
skewOMEGA = log_param.skewOMEGA./log_param.scales.time;
vreln = vecnorm(states(4:6,:) - skewOMEGA*states(1:3,:));
q = 1/2*atmo_profile.rho.*vreln.^2;
qalf = 1/2*atmo_profile.rho.*vreln.^2.*pa;
q_imp = q/6895; qalf_imp = qalf/6895;

[max_q,ind_max_q] = max(q);
[max_qalf,ind_max_qalf] = max(qalf);
ind_stg_sep = find(stages==1.5,1);
[~,ind2nd_plf_jet] = min(diff(mass(stg2_ind)));
ind_plf_jet = ind2nd_plf_jet + stg2_ind(1);

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

% velocity in surface centered local vertical local horizontal frame
vlv = dot(states(1:3,:),states(4:6,:))./rn;
vlh = sqrt(vecnorm(states(4:6,:)).^2-vlv.^2);

% mass flow rate
stg1_ind = find(stages==1); stg1p5_ind = find(stages==1.5); stg2_ind = find(stages==2);
mrate = zeros(size(states(7,:)));
mrate(stg1_ind) = min(0,[diff(states(7,stg1_ind)) 0])./dt;
mrate(stg2_ind) = min(0,[diff(states(7,stg2_ind)) 0])./dt;
patch_ind = find(t==251.02):find(t==265.02);
mrate(patch_ind) = mrate(patch_ind(1));

%% final report graphics
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % POSITION
figure; subplot_size = [7,2]; subplot_ind = 1;
subplot(subplot_size(1),subplot_size(2),subplot_ind);
plotter3d(log_x,log_param); hold on; grid off
plot3(rorb_rot(1,:),rorb_rot(2,:),rorb_rot(3,:),"Color",[0.5 0.5 0.5],"LineWidth",0.8,"LineStyle","--"); axis equal; axis off
lgd1 = legend("Trajectory","", "","Inserted Orbit","interpreter","latex","Location","southeast");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(dr/0.3048,alti/0.3048,"k","LineWidth",2); hold on; grid on
scatter(dr(ind_max_q)/0.3048,alti(ind_max_q)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(dr(ind_max_qalf)/0.3048,alti(ind_max_qalf)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","g");
scatter(dr(ind_stg_sep)/0.3048,alti(ind_stg_sep)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Downrange (ft)"); ylabel("Altitude (ft)"); axis equal
lgd1 = legend("Trajectory","Max Q", "Max Q$\alpha$","Stage Separation","interpreter","latex","Location","southeast");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,alti/0.3048,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),alti(ind_max_q)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),alti(ind_max_qalf)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","g");
scatter(t(ind_stg_sep),alti(ind_stg_sep)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Attitude (ft)");
lgd1 = legend("Trajectory","Max Q", "Max Q$\alpha$","Stage Separation","interpreter","latex","Location","southeast");
% fontsize(lgd1,20,'points'); set(lgd1,'Box','on','Color',[0.2,0.2,0.2]);
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,dr/0.3048/5280,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),dr(ind_max_q)/0.3048/5280,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),dr(ind_max_qalf)/0.3048/5280,"filled","LineWidth",2,"MarkerFaceColor","g");
scatter(t(ind_stg_sep),dr(ind_stg_sep)/0.3048/5280,"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Downrange (mi)");
lgd1 = legend("Trajectory","Max Q", "Max Q$\alpha$","Stage Separation","interpreter","latex","Location","southeast");
subplot_ind = subplot_ind+1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % END POSITION

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % VELOCITY
subplot(subplot_size(1),subplot_size(2),subplot_ind);
plotter3d(log_x,log_param); hold on; grid off; axis off
legend("Trajectory","Velocity Vector Visualization","interpreter","latex","Location","southwest");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,pa/pi*180,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),pa(ind_max_q)/pi*180,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),pa(ind_max_qalf)/pi*180,"filled","LineWidth",2,"MarkerFaceColor","g");
scatter(t(ind_stg_sep),pa(ind_stg_sep)/pi*180,"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Flight Path Angle (deg)");
legend("Trajectory","Max Q", "Max Q$\alpha$","Stage Separation","interpreter","latex","Location","northeast");
subplot_ind = subplot_ind+1;

% subplot(subplot_size(1),subplot_size(2),subplot_ind);
% plot(t,vn/0.3048,"k","LineWidth",2); hold on; grid on
% scatter(t(ind_max_q),vn(ind_max_q)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","b");
% scatter(t(ind_max_qalf),vn(ind_max_qalf)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","g");
% scatter(t(ind_stg_sep),vn(ind_stg_sep)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","r");
% xlabel("Time since Release (s)"); ylabel("Speed (ft/s)");
% legend("Trajectory","Max Q", "Max Q$\alpha$","Stage Separation","interpreter","latex","Location","southeast");
% subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,vlv/0.3048,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),vlv(ind_max_q)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),vlv(ind_max_qalf)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","g");
scatter(t(ind_stg_sep),vlv(ind_stg_sep)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Climb Rate (ft/s)");
legend("Trajectory","Max Q", "Max Q$\alpha$","Stage Separation","interpreter","latex","Location","northeast");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,vlh/0.3048,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),vlh(ind_max_q)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),vlh(ind_max_qalf)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","g");
scatter(t(ind_stg_sep),vlh(ind_stg_sep)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Downrange Speed (ft/s)");
legend("Trajectory","Max Q", "Max Q$\alpha$","Stage Separation","interpreter","latex","Location","southeast");
subplot_ind = subplot_ind+1;

% subplot(subplot_size(1),subplot_size(2),subplot_ind);
% plot(t,states(4,:)/0.3048,"k","LineWidth",2); hold on; grid on
% plot(t,states(5,:)/0.3048,"k--","LineWidth",2);
% plot(t,states(6,:)/0.3048,"k:","LineWidth",2);
% xlabel("Time since Release (s)"); ylabel("Velocity Component in ECI (ft/s)");
% legend("$v_x$","$v_y$", "$v_z$","interpreter","latex","Location","southeast");
% subplot_ind = subplot_ind+1;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % END VELOCITY

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % MASS
subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,states(7,:)*2.20462262,"k","LineWidth",2); hold on; grid on
scatter(t(ind_stg_sep),states(7,ind_stg_sep)*2.20462262,"filled","LineWidth",2,"MarkerFaceColor","r");
scatter(t(ind_plf_jet),states(7,ind_plf_jet)*2.20462262,"filled","LineWidth",2,"MarkerFaceColor","m");
xlabel("Time since Release (s)"); ylabel("Mass (lb)");
legend("Mass","Stage Separation","Payload Fairing Separation","interpreter","latex","Location","northeast");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,mrate*2.20462262,"k","LineWidth",2); hold on; grid on
scatter(t(ind_stg_sep),mrate(ind_stg_sep)*2.20462262,"filled","LineWidth",2,"MarkerFaceColor","r");
scatter(t(ind_plf_jet),mrate(ind_plf_jet)*2.20462262,"filled","LineWidth",2,"MarkerFaceColor","m");
xlabel("Time since Release (s)"); ylabel("Mass Flow Rate (lb/s)");
legend("Mass Flow Rate","Stage Separation","Payload Fairing Separation","interpreter","latex","Location","southeast");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,rcg(1,:)/0.3048,"k","LineWidth",2); hold on; grid on
scatter(t(ind_stg_sep),rcg(1,ind_stg_sep)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","r");
scatter(t(ind_plf_jet),rcg(1,ind_plf_jet)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","m");
xlabel("Time since Release (s)"); ylabel("Center of Mass Location (ft)"); ylim([-param.vehicle_sizing.lengths(3) 0]);
legend("$r_{cg}$","Stage Separation","Payload Fairing Separation","interpreter","latex","Location","northwest");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t(stg1_ind),gs_load(stg1_ind),"k","LineWidth",2); hold on
plot(t(stg2_ind),gs_load(stg2_ind),"k","LineWidth",2); grid on
xlabel("Time since Release (s)"); ylabel("Unsigned Payload G Loading (gs)");
% legend("1st Stage Propulsion","2nd Stage Propulsion","interpreter","latex");
subplot_ind = subplot_ind+1;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % END MASS

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % MAX Q
subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,q_imp,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),q_imp(ind_max_q),"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),q_imp(ind_max_qalf),"filled","LineWidth",2,"MarkerFaceColor","r");
legend("Max Q Time History","Max Q", "Max Q$\alpha$","interpreter","latex");
ylabel("Dynamics Pressure (psi)");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,qalf_imp,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),qalf_imp(ind_max_q),"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),qalf_imp(ind_max_qalf),"filled","LineWidth",2,"MarkerFaceColor","r");
legend("Max Q$\alpha$ Time History","Max Q", "Max Q$\alpha$","interpreter","latex");
ylabel("$Q\alpha$ Pressure (psi)"); xlabel("Time since Release (s)");
subplot_ind = subplot_ind+1;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % END MAX Q

%% state plotter

figure; subplot_size = [2,3]; subplot_ind = 1;
subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,alti/1000,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),alti(ind_max_q)/1000,"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),alti(ind_max_qalf)/1000,"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Vehicle Attitude (km)");
legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,vn,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),vn(ind_max_q),"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),vn(ind_max_qalf),"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,dr/1000,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Downrange (km)");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,states(7,:),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,pa*180/pi,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),pa(ind_max_q)*180/pi,"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),pa(ind_max_qalf)*180/pi,"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Pitch Angle (deg)");
legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");
subplot_ind = subplot_ind+1;

figure; subplot_size = [2,3]; subplot_ind = 1;
subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,alti/0.3048,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),alti(ind_max_q)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),alti(ind_max_qalf)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","g");
xlabel("Time since Release (s)"); ylabel("Vehicle Attitude (ft)");
lgd1 = legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");
% fontsize(lgd1,20,'points'); set(lgd1,'Box','on','Color',[0.2,0.2,0.2]);
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(dr/0.3048,alti/0.3048,"k","LineWidth",2); grid on
xlabel("Vehicle Downrange (ft)"); ylabel("Vehicle Altitude (ft)"); axis equal
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,dr/0.3048/5280,"k","LineWidth",2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Downrange (mi)");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,vn/0.3048,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),vn(ind_max_q)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),vn(ind_max_qalf)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","g");
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (ft/s)");
lgd2 = legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");
% fontsize(lgd2,20,'points'); set(lgd2,'Box','on','Color',[0.2,0.2,0.2]);
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,states(7,:)*2.20462262,"k","LineWidth",2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (lb)");
subplot_ind = subplot_ind+1;

subplot(subplot_size(1),subplot_size(2),subplot_ind);
plot(t,pa*180/pi,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),pa(ind_max_q)*180/pi,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),pa(ind_max_qalf)*180/pi,"filled","LineWidth",2,"MarkerFaceColor","g");
xlabel("Time since Release (s)"); ylabel("Pitch Angle (deg)");
lgd4 = legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");
% fontsize(lgd4,20,'points'); set(lgd4,'Box','on','Color',[0.2,0.2,0.2]); plot_darkmode
subplot_ind = subplot_ind+1;

%% control plotter
figure;
subplot(2,2,1); plot(t,tvc(1,:),"r","LineWidth",1.2); hold on; grid on
plot(t,tvc(2,:),"b","LineWidth",1.2); plot(t,tvc(3,:),"k","LineWidth",1.2);
xlabel("Time since Release (s)"); ylabel("Thrust Vector Control");
subplot(2,2,2); plot(t,throtl*100,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Throttle (%)");
subplot(2,2,4); plot(t,tvca*180/pi,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("TVC Body Angle (deg)");

%% g acceleration/loading plotter
figure; subplot(2,1,1);
plot(t(stages==1),gs_acce(stages==1),"b","LineWidth",1.2); hold on
plot(t(stages==2),gs_acce(stages==2),"r","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Unsigned G Acceleration (gs)");
legend("1st Stage Burn","2nd Stage Burn","interpreter","latex");
subplot(2,1,2);
plot(t(stages==1),gs_load(stages==1),"b","LineWidth",1.2); hold on
plot(t(stages==2),gs_load(stages==2),"r","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Unsigned G Loading (gs)");
legend("1st Stage Burn","2nd Stage Burn","interpreter","latex");

%% flight limit load factor
fllf = [0 0.5 0.5 1.5 1.5 1 1 0; 6 6 4 3 -1.5 -1.5 -2 -2];
fllf_falcon9 = [0 0.5 0.5 2 2 0.5 0.5 0; 6 6 4 3.5 -1.5 -1.5 -2 -2];
fllf = [fllf flip([-fllf(1,:);fllf(2,:)],2)];
fllf_falcon9 = [fllf_falcon9 flip([-fllf_falcon9(1,:);fllf_falcon9(2,:)],2)];
figure;
scatter(gacce_lateral(stages==1),gacce_axial(stages==1),"b","filled"); hold on
scatter(gacce_lateral(stages==2),gacce_axial(stages==2),"r","filled"); grid on
plot(fllf_falcon9(1,:),fllf_falcon9(2,:),"Color",[0.5 0.5 0.5],"LineWidth",1.5,"LineStyle","--");
plot(fllf(1,:),fllf(2,:),"Color",[1 1 1],"LineWidth",2);
xlim([-2.5 2.5]); ylim([-2.5 6.5]);
xlabel("Lateral Acceleration (g)"); ylabel("Axial Acceleration (g)");
lgd = legend("1st Stage Burn","2nd Stage Burn","Falcon 9","Thumper","interpreter","latex");
fontsize(lgd,20,'points'); plot_darkmode; set(lgd,'Box','on','Color',[0.2,0.2,0.2]);

%% dynamic pressure plot
figure; subplot(2,1,1); plot(t,q,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),q(ind_max_q),"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),q(ind_max_qalf),"filled","LineWidth",2,"MarkerFaceColor","r");
ylabel("Dynamics Pressure (Pa)");
subplot(2,1,2); plot(t,qalf,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),qalf(ind_max_q),"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),qalf(ind_max_qalf),"filled","LineWidth",2,"MarkerFaceColor","r");
ylabel("$Q\alpha$ Pressure (Pa)"); xlabel("Time since Release (s)");

figure; subplot(2,1,1); plot(t,q_imp,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),q_imp(ind_max_q),"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),q_imp(ind_max_qalf),"filled","LineWidth",2,"MarkerFaceColor","r");
legend("Max Q Time History","Max Q", "Max Q$\alpha$","interpreter","latex");
ylabel("Dynamics Pressure (psi)");
subplot(2,1,2); plot(t,qalf_imp,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),qalf_imp(ind_max_q),"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),qalf_imp(ind_max_qalf),"filled","LineWidth",2,"MarkerFaceColor","r");
legend("Max Q$\alpha$ Time History","Max Q", "Max Q$\alpha$","interpreter","latex");
ylabel("$Q\alpha$ Pressure (psi)"); xlabel("Time since Release (s)");

figure; semilogy(t,q_imp,"k","LineWidth",1.2); hold on
semilogy(t,qalf_imp,"r","LineWidth",1.2); grid on
ylabel("Aerodynamic Loading log(psi)"); xlabel("Time since Release (s)");
legend("Q","$Q\alpha$","interpreter","latex");

save("continous_traj.mat","t","states","ctrls");
