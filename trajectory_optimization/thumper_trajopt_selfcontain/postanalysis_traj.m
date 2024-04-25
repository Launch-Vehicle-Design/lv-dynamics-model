clear; clc; close all
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

file_name = "thumper_straj_cg3dof.mat";
if exist(file_name,"file")
    load(file_name,"log_x","log_param");
end

% set mission time since the drop
t = 0:0.02:520;
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

%% state plotter
figure;
subplot(2,2,1); plot(t,alti/1000,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),alti(ind_max_q)/1000,"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),alti(ind_max_qalf)/1000,"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Vehicle Attitude (km)");
legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");

subplot(2,2,2); plot(t,vn,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),vn(ind_max_q),"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),vn(ind_max_qalf),"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");

subplot(2,2,3); plot(t,states(7,:),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");

subplot(2,2,4); plot(t,pa*180/pi,"k","LineWidth",1.2); hold on; grid on
scatter(t(ind_max_q),pa(ind_max_q)*180/pi,"filled","LineWidth",2,"MarkerFaceColor","m");
scatter(t(ind_max_qalf),pa(ind_max_qalf)*180/pi,"filled","LineWidth",2,"MarkerFaceColor","r");
xlabel("Time since Release (s)"); ylabel("Pitch Angle (deg)");
legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");

figure;
subplot(2,2,1); plot(t,alti/0.3048,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),alti(ind_max_q)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),alti(ind_max_qalf)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","g");
xlabel("Time since Release (s)"); ylabel("Vehicle Attitude (ft)");
lgd1 = legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");
% fontsize(lgd1,20,'points'); set(lgd1,'Box','on','Color',[0.2,0.2,0.2]);

subplot(2,2,2); plot(t,vn/0.3048,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),vn(ind_max_q)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),vn(ind_max_qalf)/0.3048,"filled","LineWidth",2,"MarkerFaceColor","g");
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (ft/s)");
lgd2 = legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");
% fontsize(lgd2,20,'points'); set(lgd2,'Box','on','Color',[0.2,0.2,0.2]);

subplot(2,2,3); plot(t,states(7,:)*2.20462262,"k","LineWidth",2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (lb)");

subplot(2,2,4); plot(t,pa*180/pi,"k","LineWidth",2); hold on; grid on
scatter(t(ind_max_q),pa(ind_max_q)*180/pi,"filled","LineWidth",2,"MarkerFaceColor","b");
scatter(t(ind_max_qalf),pa(ind_max_qalf)*180/pi,"filled","LineWidth",2,"MarkerFaceColor","g");
xlabel("Time since Release (s)"); ylabel("Pitch Angle (deg)");
lgd4 = legend("Trajectory","Max Q", "Max Q$\alpha$","interpreter","latex");
% fontsize(lgd4,20,'points'); set(lgd4,'Box','on','Color',[0.2,0.2,0.2]); plot_darkmode

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
