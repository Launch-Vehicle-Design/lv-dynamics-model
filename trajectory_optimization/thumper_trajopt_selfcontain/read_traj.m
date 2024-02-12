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
for i = 1:length(t)
    cont_traj = cont3d_traj_interp(t(i),log_x,log_param);
    states(:,i) = cont_traj.state; ctrls(:,i) = cont_traj.ctrl;
end

earth_radius = log_param.earthR*log_param.scales.length;
rn = vecnorm(states(1:3,:)); vn = vecnorm(states(4:6,:)); alti = rn-earth_radius;
tvc = ctrls(1:3,:); throtl = ctrls(4,:);
pa = acos(dot(states(4:6,:)./vn,tvc));

figure;
subplot(2,2,1); plot(t,alti/1000,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Attitude (km)");
subplot(2,2,2); plot(t,vn,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
subplot(2,2,3); plot(t,states(7,:),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");
subplot(2,2,4); plot(t,pa*180/pi,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Pitch Angle (deg)");

figure;
subplot(2,2,1); plot(t,tvc(1,:),"r","LineWidth",1.2); grid on
hold on; plot(t,tvc(2,:),"b","LineWidth",1.2); plot(t,tvc(3,:),"k","LineWidth",1.2);
xlabel("Time since Release (s)"); ylabel("Thrust Vector Control");
subplot(2,2,2); plot(t,throtl*100,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Throttle (%)");

atmo_profile = atmo(alti);
q = 1/2*atmo_profile.rho.*vn.^2;
qalf_uplim = 1/2*atmo_profile.rho.*vn.^2.*pa;

figure; subplot(1,2,1); plot(t,q,"k","LineWidth",1.2); grid on
ylabel("Dynamics Pressure (Pa)");
subplot(1,2,2); plot(t,qalf_uplim,"k","LineWidth",1.2); grid on
ylabel("Q $\alpha$ Pressure (Pa)"); xlabel("Time since Release (s)");

save("continous_traj.mat","t","states","ctrls");
