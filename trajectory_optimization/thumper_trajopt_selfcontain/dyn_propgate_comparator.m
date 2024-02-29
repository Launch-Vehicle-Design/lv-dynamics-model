clear; clc; close all
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

file_name = "thumper_straj_cgp3dof.mat";
if exist(file_name,"file")
    load(file_name,"log_x","log_param");
end

% set mission time since the drop
dt = 1;
tf = 700;
t = 0:dt:tf;
init_cond = log_param.init_cond;
x0 = log_param.x0;

states = nan([log_param.nstate,length(t)]);
ctrls = nan([log_param.nctrl,length(t)]);
forces = nan([9,length(t)]);
stages = nan([1,length(t)]);

for i = 1:length(t)
    cont_traj = cont3d_traj_interp(t(i),log_x,log_param);
    states(:,i) = cont_traj.state; ctrls(:,i) = cont_traj.ctrl;
    forces(:,i) = cont_traj.force; stages(i) = cont_traj.stage;
end

% dynamic propagation
dt_dyn = dt/cont_traj.scales.time;
t_dyn = t./cont_traj.scales.time;
states_dyn = nan([log_param.nstate,length(t_dyn)]);

dynfunc_1st_burn = @(t,x,u) rocket(t,x,u,log_param,1);
dynfunc_2nd_burn = @(t,x,u) rocket(t,x,u,log_param,2);

current_x = init_cond; enddrop_ind = find(stages==1,1)-1;
record = rk4sim(dynfunc_1st_burn,current_x,dt_dyn,t_dyn(enddrop_ind),ctrls(:,1:enddrop_ind));
states_dyn(:,1:enddrop_ind) = record.x.*cont_traj.scales.state;

current_x = [record.x(1:3,end); x0(4:6); record.x(7,end)]; end1st_ind = find(stages==2,1)-1;
record = rk4sim(dynfunc_1st_burn,current_x,dt_dyn,t_dyn(end1st_ind)-t_dyn(enddrop_ind+1),ctrls(:,enddrop_ind+1:end1st_ind));
states_dyn(:,enddrop_ind+1:end1st_ind) = record.x.*cont_traj.scales.state;

current_x = [record.x(1:6,end); log_param.m02];
record = rk4sim(dynfunc_2nd_burn,current_x,dt_dyn,t_dyn(end)-t_dyn(end1st_ind+1),ctrls(:,end1st_ind+1:end));
states_dyn(:,end1st_ind+1:end) = record.x.*cont_traj.scales.state;

% state plotter
earth_radius = log_param.earthR*log_param.scales.length;
rn = vecnorm(states(1:3,:)); vn = vecnorm(states(4:6,:)); alti = rn-earth_radius;
tvc = ctrls(1:3,:); throtl = ctrls(4,:); pa = acos(dot(states(4:6,:)./vn,tvc));
rn_dyn = vecnorm(states_dyn(1:3,:)); vn_dyn = vecnorm(states_dyn(4:6,:));
alti_dyn = rn_dyn-earth_radius; pa_dyn = acos(dot(states_dyn(4:6,:)./vn_dyn,tvc));

figure;
subplot(2,2,1);
plot(t,alti/1000,"k","LineWidth",1.2); hold on
plot(t,alti_dyn/1000,"b","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Attitude (km)");
subplot(2,2,2);
plot(t,vn,"k","LineWidth",1.2); hold on
plot(t,vn_dyn,"b","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
subplot(2,2,3);
plot(t,states(7,:),"k","LineWidth",1.2); hold on
plot(t,states_dyn(7,:),"b","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");
subplot(2,2,4);
plot(t,pa*180/pi,"k","LineWidth",1.2); hold on
plot(t,pa_dyn*180/pi,"b","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("TVC Angle (deg)");
