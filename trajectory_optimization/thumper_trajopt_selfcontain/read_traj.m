clear; clc; close all
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

file_name = "thumper_straj_cg3dof.mat";
if exist(file_name,"file")
    load(file_name,"log_x","log_param");
end

% set mission time since the drop
t = 0:0.1:500;
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

force_acce = forces(4:6,:)./states(7,:) + forces(7:9,:)./states(7,:);
acce = forces(1:3,:) + force_acce;
gs_acce = vecnorm(acce)/(log_param.g0*log_param.scales.acceleration);
gs_load = vecnorm(force_acce)/(log_param.g0*log_param.scales.acceleration);

% state plotter
figure;
subplot(2,2,1); plot(t,alti/1000,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Attitude (km)");
subplot(2,2,2); plot(t,vn,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
subplot(2,2,3); plot(t,states(7,:),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");
subplot(2,2,4); plot(t,pa*180/pi,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Pitch Angle (deg)");

% control plotter
figure;
subplot(2,2,1); plot(t,tvc(1,:),"r","LineWidth",1.2); hold on; grid on
plot(t,tvc(2,:),"b","LineWidth",1.2); plot(t,tvc(3,:),"k","LineWidth",1.2);
xlabel("Time since Release (s)"); ylabel("Thrust Vector Control");
subplot(2,2,2); plot(t,throtl*100,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Throttle (%)");

% g acceleration/loading plotter
figure; subplot(2,1,1);
plot(t(stages==1),gs_acce(stages==1),"b","LineWidth",1.2); hold on
plot(t(stages==2),gs_acce(stages==2),"r","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Unsigned G Acceleration (gs)");
legend("1st Stage", "2nd Stage");
subplot(2,1,2);
plot(t(stages==1),gs_load(stages==1),"b","LineWidth",1.2); hold on
plot(t(stages==2),gs_load(stages==2),"r","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Unsigned G Loading (gs)");
legend("1st Stage", "2nd Stage");

% dynamic pressure
atmo_profile = atmo(alti);
q = 1/2*atmo_profile.rho.*vn.^2;
qalf_uplim = 1/2*atmo_profile.rho.*vn.^2.*pa;

figure; subplot(2,1,1); plot(t,q,"k","LineWidth",1.2); grid on
ylabel("Dynamics Pressure (Pa)");
subplot(2,1,2); plot(t,qalf_uplim,"k","LineWidth",1.2); grid on
ylabel("Q$\alpha$ Pressure (Pa)"); xlabel("Time since Release (s)");

figure; semilogy(t,q,"k","LineWidth",1.2); hold on
semilogy(t,qalf_uplim,"r","LineWidth",1.2); grid on
ylabel("Q$\alpha$ Pressure (Pa)"); xlabel("Time since Release (s)");

save("continous_traj.mat","t","states","ctrls");
