function [u_tvc,u_gf] = mpctvc1st(curr_x,curr_u,traj,param,phase)
mpc.n_hrzn = 10;
mpc.n_state = 14; mpc.n_ctrl = 6;
mpc.ntot = mpc.n_state + mpc.n_ctrl;
mpc.state_ind = 1:14;

mpc.dt = param.mpc_dt;
mpc.t = traj.t;
mpc.phase = phase;

[~,mpc.start_ind] = min(vecnorm(traj.states(1:3,:)-curr_x(1:3)));
mpc.ref_states = traj.states(:,mpc.start_ind:mpc.start_ind+mpc.n_hrzn-1);
mpc.init_x = curr_x;
mpc.init_u = curr_u;
mpc.init_C = EP2C(curr_x(7:10));

% construct initial optimization variable and bounds
mpc.ref_traj = [mpc.ref_states(1:6,:); repmat(curr_x(7:13),[1,mpc.n_hrzn]); mpc.ref_states(7,:)];
delta_x = repmat(curr_x(mpc.state_ind),[1,mpc.n_hrzn]) - mpc.ref_traj;
init_x =  reshape([delta_x' rand([mpc.n_hrzn,mpc.n_ctrl])],[mpc.ntot*mpc.n_hrzn,1]);

% tvc actuation bounds
tvc_bounds = zeros([2*mpc.n_hrzn,2]);
tvc_bounds(:,1) = param.tvc_limit(1);
tvc_bounds(:,2) = param.tvc_limit(2);
% grid fin actuation bounds
gf_bounds = zeros([4*mpc.n_hrzn,2]);
gf_bounds(:,1) = param.gf_limit(1);
gf_bounds(:,2) = param.gf_limit(2);
% lower bounds for the optimization valuables
low_bound = [-Inf*ones(6*mpc.n_hrzn,1); -1*ones(4*mpc.n_hrzn,1); -Inf*ones(3*mpc.n_hrzn,1); tvc_bounds(:,1); gf_bounds(:,1)];
upp_bound = [Inf*ones(6*mpc.n_hrzn,1); ones(4*mpc.n_hrzn,1); Inf*ones(3*mpc.n_hrzn,1); tvc_bounds(:,2); gf_bounds(:,2)];

% linear state propagation equality constraint
% [mpc.A,mpc.B] = lv1st3dof(curr_x,curr_u,param);
[mpc.A,mpc.B] = lv1st([mpc.ref_traj(:,1);curr_x(15:16)],curr_u,param);
comb_mat = [mpc.A mpc.B; zeros(mpc.n_ctrl,mpc.ntot)];
comb_mat_d = expm(comb_mat*mpc.dt);
mpc.Ad = comb_mat_d(1:mpc.n_state,1:mpc.n_state);
mpc.Bd = comb_mat_d(1:mpc.n_state,mpc.n_state+1:end);

sparse_a_mat = [zeros(1,mpc.n_hrzn); [eye(mpc.n_hrzn-1) zeros(mpc.n_hrzn-1,1)]];
Aeq = [-eye(mpc.n_hrzn*mpc.n_state)+kron(sparse_a_mat,mpc.Ad) kron(eye(mpc.n_hrzn),mpc.Bd)];
Beq = [-mpc.Ad*delta_x(:,1); zeros((mpc.n_hrzn-1)*mpc.n_state,1)];

% setup cost and constraint functions
cost_func = @(x) cost(x,mpc,traj,param);
% initial state cost and constraint function
init_cost = cost_func(init_x);

options = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point',...
    "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',1e8,'MaxIterations',100,'EnableFeasibilityMode',true);
x = fmincon(cost_func,init_x,[],[],Aeq,Beq,low_bound,upp_bound,[],options);
u = x(mpc.n_hrzn*mpc.n_state+1:mpc.n_hrzn:mpc.ntot*mpc.n_hrzn);

    function val = cost(x,mpc,traj,param)
        persistent func_count;
        if isempty(func_count)
            func_count = 0;
        end

        n = mpc.n_hrzn; ntot = mpc.ntot;
        n_ctrl = mpc.n_ctrl; n_state = mpc.n_state;
        states = reshape(x(1:n*ntot),[n,ntot])';
        val = 0.01*sum(vecnorm(states(1:3,:))) + sum(vecnorm(states(4:6,:))) + ...
            sum(vecnorm(states(7:10,:))) + sum(vecnorm(states(11:13,:))) + ...
            0.1*sum(abs(states(14,:)));

        if mod(func_count,5000) == 0
            states = reshape(x(1:n*n_state),[n,n_state])'+mpc.ref_traj;
            plotter3d(mpc.start_ind,mpc.init_x,traj,param,true,states);
            drawnow;
        end
        func_count = func_count + 1;
    end
end
