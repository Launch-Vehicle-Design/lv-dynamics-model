clear; clc; close all
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

param = sysparam()
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
cgs = nan([1,length(t)]);

for i = 1:length(t)
    cont_traj = cont3d_traj_interp(t(i),log_x,log_param);
    states(:,i) = cont_traj.state; ctrls(:,i) = cont_traj.ctrl;
    forces(:,i) = cont_traj.force; stages(i) = cont_traj.stage;
    iner = inertia(cont_traj.state(7),0,param);
    cgs(i) = iner.cg;
end

figure; plot(t,-cgs,"k","LineWidth",1.2); grid on;
xlabel("time since release, s"); ylabel("center of gravity location in nosecone frame, m");
