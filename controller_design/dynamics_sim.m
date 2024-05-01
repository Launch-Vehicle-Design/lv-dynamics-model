clear; clc; close all

addpath("rotation_toolbox\");
trajectory_file = "thumper_straj_cg3dof.mat";
interp_traj_file = "interpolated_traj.mat";
featured_dt = 0.01; featured_tf = 510;
t = 0:featured_dt:featured_tf;

load(trajectory_file,"log_x","log_param");
reinterp = ~exist(interp_traj_file,"file");
if ~reinterp
    load(interp_traj_file);
    if traj.dt ~= featured_dt || traj.tf ~= featured_tf
        reinterp = true;
    end
end
if reinterp
    % read and interpolate the trajectory
    states = nan([log_param.nstate,length(t)]);
    ctrls = nan([log_param.nctrl,length(t)]);
    forces = nan([9,length(t)]);
    stages = nan([1,length(t)]);
    for i = 1:length(t)
        cont_traj = cont3d_traj_interp(t(i),log_x,log_param);
        states(:,i) = cont_traj.state; ctrls(:,i) = cont_traj.ctrl;
        forces(:,i) = cont_traj.force; stages(i) = cont_traj.stage;
    end
    traj.dt = featured_dt; traj.tf = featured_tf; traj.t = t; traj.states = states;
    traj.ctrls = ctrls; traj.forces = forces; traj.stages = stages;
    save(interp_traj_file,"traj");
end

% define system parameters
param = sysparam();
param.mpc_dt = featured_dt;

% extract the initial condition from logged parameter
state_scaling = [log_param.scales.length*ones([3,1]);
    log_param.scales.speed*ones([3,1]); log_param.scales.mass];
init_cond_3dof = log_param.init_cond.*state_scaling;

% re-state the initial condition for 6dof
init_r = init_cond_3dof(1:3); init_v = init_cond_3dof(4:6);
e_b1 = init_v/norm(init_v);
e_b2 = cross(e_b1,init_r)/norm(cross(e_b1,init_r));
e_b3 = cross(e_b1,e_b2);
init_C = [e_b1 e_b2 e_b3]; init_q = C2EP(init_C);
init_omega = zeros([3,1]);
init_m = init_cond_3dof(7); init_mox = param.m.ox2nd;
init_mfuel = param.m.fuel2nd;
init_oxang = zeros([2,1]); init_fuelang = zeros([2,1]);
init_doxang = zeros([2,1]); init_dfuelang = zeros([2,1]);
init_x = [init_r; init_v; init_q; init_omega;
    init_m; init_mox; init_mfuel; init_oxang;
    init_fuelang; init_doxang; init_dfuelang];

% simulation
mpc_hrzn_len = 20;
curr_x = init_x;
u_acs = zeros([8,1]); u_gf = zeros([4,1]);
for i = 1:size(t,2)
    curr_t = t(i);
    curr_C = EP2C(curr_x(7:10));
    if traj.stages(i) == 0
        u = [zeros([3,1]); u_acs; u_gf];
    elseif traj.stages(i) == 1
        [~,start_ind] = min(vecnorm(curr_x(1:3)-traj.states(1:3,:)));
        traj.ref_ind = start_ind:start_ind+mpc_hrzn_len-1;
        traj.ref_ind(traj.ref_ind>size(t,2)) = size(t,2);
        % extract from trajectory profile - naive control
        tvc_iner = traj.ctrls(1:3,i);
        tvc_body = curr_C'*tvc_iner;
        p = -asin(tvc_body(3));
        y = atan(tvc_body(2)/tvc_body(1));
        curr_u = [y; max(min(p,param.tvc_limit(2)),param.tvc_limit(1)); u_gf];

        % 1st stage main engine TVC controller
        [u_tvc,u_gf] = mpctvc1st(curr_x,curr_u,traj,param,traj.stages(i),mpc_hrzn_len);
        u = [u_tvc; traj.ctrls(4,i); u_acs; u_gf];
    end

    dynFunc = @(t,x,u) rocket6dof(t,x,u,param,traj.stages(i));
    plotter3d(i,curr_x,traj,param); drawnow;

    x_hist(:,i) = curr_x; u_hist(:,i) = u;
    curr_x = singleRK4(dynFunc,curr_t,featured_dt,curr_x,u);

    % result = ode45(dynFunc,[0 featured_dt],curr_x);
    % curr_x = result.y(:,end);
end
