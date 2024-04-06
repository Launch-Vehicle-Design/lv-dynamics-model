function u = nmpc(curr_t,curr_x,traj,param,phase)
mpc.n_hrzn = 10;
mpc.n_state = 24; mpc.n_ctrl = 11;
mpc.ntot = mpc.n_state + mpc.n_ctrl;
mpc.dt = param.mpc_dt;
mpc.t = traj.t;
mpc.phase = phase;

% scales
mpc.scales.length = param.earthR;
mpc.scales.speed = sqrt(param.mu/param.earthR);
mpc.scales.time = mpc.scales.length/mpc.scales.speed;
mpc.scales.acce = mpc.scales.speed/mpc.scales.time;
mpc.scales.mass = param.m0;
mpc.scales.force = mpc.scales.mass*mpc.scales.acce;
mpc.scales.traj_state = [mpc.scales.length*ones([3,1]);
    mpc.scales.speed*ones([3,1]); mpc.scales.mass];
mpc.scales.opt_state = [mpc.scales.length*ones([3,1]);
    mpc.scales.speed*ones([3,1]); ones([7,1]);
    mpc.scales.mass*ones([3,1]); ones(8,1)];

mpc.start_ind = find(traj.t == curr_t);
mpc.ref_states = traj.states(:,mpc.start_ind:mpc.start_ind+mpc.n_hrzn-1);
mpc.init_x = curr_x;

% pre-propagate the target location in MCO (ground truth)
% tspan = (1:mpc.n_hrzn)*mpc.dt;
% rocket_ode = @(t,x) rocket6dof(t,x,zeros([mpc.n_ctrl,1]),param,phase);
% [~,mpc.init_rollout] = ode45(rocket_ode,tspan,curr_x);

% construct initial optimization variable and bounds
init_x = reshape([repmat((curr_x./mpc.scales.opt_state)',[mpc.n_hrzn,1]) zeros([mpc.n_hrzn,mpc.n_ctrl])],[mpc.ntot*mpc.n_hrzn,1]);

% tvc actuation bounds
tvc_bounds = zeros([2*mpc.n_hrzn,2]);
thr_bounds = zeros([mpc.n_hrzn,2]);
if phase == 1
    tvc_bounds(:,1) = param.tvc_limit(1);
    tvc_bounds(:,2) = param.tvc_limit(2);
    thr_bounds(:,1) = param.T_limit(1);
    thr_bounds(:,2) = param.T_limit(2);
elseif phase == 2 || phase == 3
    tvc_bounds(:,1) = param.litvc_limit(1);
    tvc_bounds(:,2) = param.litvc_limit(2);
    thr_bounds(:,1) = param.T_limit(3);
    thr_bounds(:,2) = param.T_limit(4);
end
% lower bounds for the optimization valuables
low_bound = [-1.2*repmat(vecnorm(mpc.ref_states(1:3,:)),[1,3])'/mpc.scales.length;
    -1.2*repmat(vecnorm(mpc.ref_states(4:6,:)),[1,3])'/mpc.scales.speed;
    -1*ones([4*mpc.n_hrzn,1]); -Inf*ones([3*mpc.n_hrzn,1]);
    zeros([3*mpc.n_hrzn,1]);
    -pi*ones([4*mpc.n_hrzn,1]); -Inf*ones([4*mpc.n_hrzn,1]);
    tvc_bounds(:,1); thr_bounds(:,1);
    -1*ones([8*mpc.n_hrzn,1])];
upp_bound = [1.2*repmat(vecnorm(mpc.ref_states(1:3,:)),[1,3])'/mpc.scales.length;
    1.2*repmat(vecnorm(mpc.ref_states(4:6,:)),[1,3])'/mpc.scales.speed;
    ones([4*mpc.n_hrzn,1]); Inf*ones([3*mpc.n_hrzn,1]);
    param.m0*ones([mpc.n_hrzn,1])/mpc.scales.mass;
    param.m.ox2nd*ones([mpc.n_hrzn,1])/mpc.scales.mass;
    param.m.fuel2nd*ones([mpc.n_hrzn,1])/mpc.scales.mass;
    pi*ones([4*mpc.n_hrzn,1]); Inf*ones([4*mpc.n_hrzn,1]);
    tvc_bounds(:,2); thr_bounds(:,2);
    ones([8*mpc.n_hrzn,1])];

% setup cost and constraint functions
cost_func = @(x) cost(x,mpc);
nonc_func = @(x) nonl_constraint(x,mpc,param);

% initial state cost and constraint function
init_cost = cost_func(init_x);
[init_c, init_ceq] = nonc_func(init_x);

options = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point',...
    "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',1e8,'MaxIterations',1000,'EnableFeasibilityMode',true);
x = fmincon(cost_func,init_x,[],[],[],[],low_bound,upp_bound,nonc_func,options);
u = x(mpc.n_hrzn*mpc.n_state+(1:mpc.n_hrzn:mpc.n_hrzn*mpc.n_ctrl));

    function val = cost(x,mpc)
        % cost weighting
        w_dr = 0.01; w_fr = 1; w_dv = 1; w_fv = 1; w_m = 1e2;
        
        ref_states = mpc.ref_states./mpc.scales.traj_state;
        states = reshape(x(1:mpc.n_state*mpc.n_hrzn),[mpc.n_hrzn,mpc.n_state])';
        val = w_dr*sum(vecnorm(states(1:3,:) - ref_states(1:3,:))) + ...
            w_dv*sum(vecnorm(states(4:6,:) - ref_states(4:6,:))) + ...
            w_fr*vecnorm(states(1:3,end) - ref_states(1:3,end)) + ...
            w_fv*vecnorm(states(4:6,end) - ref_states(4:6,end)) - ...
            w_m*states(14,end);
    end

    function [c,ceq] = nonl_constraint(x,mpc,param)
        persistent func_count;
        if isempty(func_count) 
            func_count = 0;
        end

        % dynamics feasibility by Euler forward rollout
        n = mpc.n_hrzn; ntot = mpc.ntot;
        n_state = mpc.n_state; n_ctrl = mpc.n_ctrl;
        initx = mpc.init_x./mpc.scales.opt_state;
        dt = mpc.dt./mpc.scales.time;
        states = [initx reshape(x(1:n_state*n),[n,n_state])'];
        % condition the control input
        contrl_signal = x(n*n_state+1:n*ntot);
        contrl = reshape(contrl_signal,[n,n_ctrl])';

        dyn_states = states(:,1:end-1).*mpc.scales.opt_state;
        dyn_dxdt = rocket6dof([],dyn_states,contrl,param,mpc.phase*ones([1,n]));
        dxdt = dyn_dxdt.*(mpc.scales.opt_state/mpc.scales.time);
        
        ref_states = mpc.ref_states./mpc.scales.traj_state;
        ceq = [reshape(dxdt-1/dt*(states(:,2:end)-states(:,1:end-1)),[n_state*n,1]);
            states(1:3,end)-ref_states(1:3,end);
            states(4:6,end)-ref_states(4:6,end)];

        % % dynamics feasibility by collocation method
        % n = mpc.n_hrzn; ntot = mpc.ntot;
        % n_state = mpc.n_state; n_ctrl = mpc.n_ctrl;
        % initx = mpc.init_x./mpc.scales.opt_state;
        % dt = mpc.dt./mpc.scales.time;
        % states = [initx reshape(x(1:n_state*n),[n,n_state])'];
        % % condition the control input
        % contrl_signal = x(n*n_state+1:n*ntot);
        % contrl = [reshape(contrl_signal,[n,n_ctrl])' zeros([n_ctrl,1])];
        % 
        % dyn_states = states.*mpc.scales.opt_state;
        % dyn_dxdt = rocket6dof([],dyn_states,contrl,param,mpc.phase*ones([1,n+1]));
        % dxdt = dyn_dxdt.*(mpc.scales.opt_state/mpc.scales.time);
        % 
        % xcol = 1/2*(states(:,1:end-1)+states(:,2:end))-dt/8*diff(dxdt,1,2);
        % xcoldot = 3/2/dt*diff(states,1,2)-1/4*(dxdt(:,1:end-1)+dxdt(:,2:end));
        % ucol = 1/2*(contrl(:,1:end-1)+contrl(:,2:end));
        % 
        % dyn_xcol = xcol.*mpc.scales.opt_state;
        % dyn_dxdt_col = rocket6dof([],dyn_xcol,ucol,param,mpc.phase*ones([1,n]));
        % dxdt_col = dyn_dxdt_col.*(mpc.scales.opt_state/mpc.scales.time);
        % 
        % ref_states = mpc.ref_states./mpc.scales.traj_state;
        % ceq = [reshape(xcoldot-dxdt_col,[n_state*n,1]);
        %     states(1:3,end)-ref_states(1:3,end);
        %     states(4:6,end)-ref_states(4:6,end)];
        
        c = [];

        func_count = func_count + 1;
    end
end
