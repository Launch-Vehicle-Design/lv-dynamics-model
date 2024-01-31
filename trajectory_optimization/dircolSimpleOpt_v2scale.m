clear; clc; close all
addpath("..\parameter_functions")

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

%% Simplified Trajectory Optimization
% % Assumptions
% zero angle of attack
% zero lift generated
% 3DOF point mass launch vehicle

user_gues = true;
load_prev = true;

scales.length       = 6357e3;
scales.speed        = sqrt(3.986004418e14/6357e3);
scales.time         = scales.length/scales.speed;
scales.acceleration = scales.speed/scales.time;
scales.mass         = 2241.1532;
scales.force        = scales.mass*scales.acceleration;
scales.area         = scales.length^2;
scales.volume       = scales.area.*scales.length;
scales.density      = scales.mass/scales.volume;
scales.pressure     = scales.force/scales.area;
scales.gravparam    = scales.acceleration*scales.length^2;
param.scales = scales;

release_alti = 40000*0.3048/param.scales.length;
release_velo = 250.786/param.scales.speed;

% general parameter
param.nstate = 4;
param.nctrl = 2;

% environmental parameter
param.mu = 3.986004418e14/param.scales.gravparam;
param.earthR = 6357e3/param.scales.length;
param.h0 = 7640/param.scales.length;
param.g = @(h) 9.80665./(1+h).^2/param.scales.acceleration;
param.rho = @(h) 1.225.*exp(-h/param.h0)/param.scales.density;
param.P = @(h) 101325.*exp(-h/param.h0)/param.scales.pressure;
param.gamma = 1.4;

param.S = pi*(24*0.0254)^2/4/param.scales.area;
param.mPL = 40/param.scales.mass;
% 1st stage
param.m0 = 2241.1532/param.scales.mass;
param.ms1 = 125.2035/param.scales.mass;
param.mp1 = 1801.0036/param.scales.mass;
param.Isp1 = 276/param.scales.time;
% 2nd stage
param.m02 = 314.9462/param.scales.mass;
param.ms2 = 29.1843/param.scales.mass;
param.mp2 = 248.7618/param.scales.mass;
param.Isp2 = 369.5/param.scales.time;
% payload fairing
param.mplf = 4/param.scales.mass;

param.TtoW_1st = 1.6; param.TtoW_2nd = 0.8;
param.maxT_1st = param.TtoW_1st*param.m0*param.g(0);
param.maxT_2nd = param.TtoW_2nd*param.m02*param.g(0);

init_cond = [release_velo 0 release_alti param.m0];
% PHASE 1 - uncontrolled drop phase
drop_dt = 0.01/param.scales.time; drop_time = 5/param.scales.time;
drop_dynfun = @(t,x,u) rocket(t,x,u,param,0);
drop_record = rk4sim(drop_dynfun,init_cond,drop_dt,drop_time,[0 0]);

%% Trajectory Optimization - direct trajectory optimization
param.v_final_ref = 7754.1/param.scales.speed;
param.alt_final_ref = 251460/param.scales.length;
param.fpa_final_ref = 0;
param.plf_dropalti_ref = 140e3/param.scales.length;
% xopt represents state and control trajectory of the system
% with total length of 7*N+1
% 1:4N - state from x11...xN1:x51...x5N
% 4N+1:6N - control from u10...u1N-1,u20...u2N-1
% 6N+1:6N+2 - the step length for both 1st and 2nd stage
% 6N+3 - the coast time at stage separation

% index helpers
ntot = param.nstate+param.nctrl;

% cost function weighting
weights.mp_weight = 1e2;
weights.fpaf_weight = 10;
weights.vf_weight = 0.1;
weights.altf_weight = 0.1;
weights.dutvc_weight = 100;
weights.duT_weight = 10;

% design optimal variable bounds
tvc_bounds = [-15 15]/180*pi;
throtl_bounds = [0.8 1 0.6 1];
h1st_bounds = [0.01 5]/param.scales.time;
h2nd_bounds = [0.01 5]/param.scales.time;
stg_sep_t_bounds = [0.01 60]/param.scales.time;

N_1st = 101; N_2nd = 101;
N = N_1st+N_2nd; x0 = drop_record.x(end,:)';
x0(2) = 80*pi/180;

% initial guess propagation
param.dynfunc_1st_burn = @(t,x,u) rocket(t,x,u,param,1);
param.dynfunc_2nd_burn = @(t,x,u) rocket(t,x,u,param,2);

%% 1st initial guess - variable range with unified randomness
init_x = zeros([ntot*N+3,1]); 
param.init_ind_ptr = 0:N:ntot*N;
init_x(ntot*N+1) = 2/param.scales.time; 
init_x(ntot*N+2) = 2/param.scales.time;
init_x(ntot*N+3) = 6.7/param.scales.time;
init_x(1:param.nstate*N) = reshape(repmat([param.v_final_ref param.fpa_final_ref param.alt_final_ref param.m02-param.mp2],[N,1]),[param.nstate*N,1]);
% init_x(1:param.nstate*N) = reshape([ ...
%     (param.v_final_ref-x0(1))*rand([N,1])+x0(1) ...
%     (param.fpa_final_ref-x0(2))*rand([N,1])+x0(2) ...
%     (param.alt_final_ref-x0(3))*rand([N,1])+x0(3) ...
%     [(param.m0-param.ms1)*rand([N_1st,1])+param.ms1; ...
%     (param.m02-param.mp2)*rand([N_2nd,1])+param.ms2]
%     ],[param.nstate*N,1]);
init_x(param.nstate*N+1:(param.nstate+1)*N) = (tvc_bounds(2)-tvc_bounds(1))*rand([N,1])+tvc_bounds(1);
init_x((param.nstate+1)*N+1:ntot*N) = (throtl_bounds(2)-throtl_bounds(1))*rand([N,1])+throtl_bounds(1);

%% 2nd initial guess - non-optimal educated guess
if user_gues
current_x = x0;
for i = 1:N_1st
    init_u = [max(min(-10*(current_x(2)-80/180*pi),tvc_bounds(2)*0.99),tvc_bounds(1)*0.99) 0.8];
    dynfunc_1st = @(t,x,u) rocket(t,x,u,param,1);
    if init_x(ntot*N+1) >= 0.5/param.scales.time
        dt_1st = 0.01/param.scales.time;
        u = repmat(init_u,[ceil(init_x(ntot*N+1)/dt_1st),1]);
        record = rk4sim(dynfunc_1st,current_x',dt_1st,init_x(ntot*N+1),u);
        current_x = record.x(end,:)';
    else
        u = init_u;
        if param.m0 - current_x(4) >= param.mp1 u(2) = 0; end
        current_x = singleRK4(dynfunc_1st,0,init_x(ntot*N+1),current_x,u);
    end
    init_x(param.init_ind_ptr(1:param.nstate)+i) = current_x;
    init_x(param.init_ind_ptr(param.nstate+1:ntot)+i) = init_u;
end
dynfunc_2nd = @(t,x,u) rocket(t,x,u,param,2);
stg_sep_record = rk4sim(dynfunc_2nd,current_x',0.1,init_x(ntot*N+3),[0 0]);
current_x = stg_sep_record.x(end,:)';
init_u = [0 1];
current_x(param.nstate) = param.m02+param.mplf;
init_x(param.init_ind_ptr(1:param.nstate)+N_1st+1) = current_x;
init_x(param.init_ind_ptr(param.nstate+1:ntot)+N_1st+1) = init_u;
for i = N_1st+2:N_1st+N_2nd
    init_u = [max(min(-1*(current_x(2)-30/180*pi),tvc_bounds(2)*0.99),tvc_bounds(1)*0.99) 1];
    if init_x(ntot*N+2) >= 0.5/param.scales.time
        dt_2nd = 0.01/param.scales.time;
        u = repmat(init_u,[ceil(init_x(ntot*N+2)/dt_2nd),1]);
        record = rk4sim(dynfunc_2nd,current_x',dt_2nd,init_x(ntot*N+2),u);
        current_x = record.x(end,:)';
    else
        u = init_u;
        if param.m02 - current_x(4) >= param.mp2 u(2) = 0; end
        current_x = singleRK4(dynfunc_2nd,0,init_x(ntot*N+2),current_x,u);
    end
    init_x(param.init_ind_ptr(1:param.nstate)+i) = current_x;
    init_x(param.init_ind_ptr(param.nstate+1:ntot)+i) = init_u;
end
end

%% 3rd initial guess - load from a file with intermediate result with interpolation
file_name = "dircol_sqp_final1_scaled.mat";
if load_prev && exist(file_name,"file")
    load(file_name,"log_x","log_param");
    state_scaling = [1; 1; 1; 1]; time_scaling = 1;
    if ~contains(file_name,"scaled")
        state_scaling = [param.scales.speed; 1; param.scales.length; param.scales.mass];
        time_scaling = param.scales.time;
    end
    size_logx = size(log_x);
    log_N1st = log_param.N_1st; log_N2nd = log_param.N_2nd;
    log_N = log_N1st + log_N2nd; log_ntot = log_param.nstate + log_param.nctrl;
    log_h1st = log_x(log_ntot*log_N+1)/time_scaling;
    log_h2nd = log_x(log_ntot*log_N+2)/time_scaling;
    log_ind_ptr = log_param.init_ind_ptr;
    log_state = [x0 reshape(log_x(1:log_param.nstate*log_N),[log_N,log_param.nstate])'./state_scaling];
    log_control = [[0;1] reshape(log_x(log_param.nstate*log_N+1:log_ntot*log_N),[log_N,log_param.nctrl])'];
    % construction of initial guess
    init_x = nan([ntot*N+3,1]);
    % Hermite spline fit
    % time step size
    h_1st = log_N1st*log_h1st/N_1st;
    h_2nd = (log_N2nd-1)*log_h2nd/N_2nd;
    % derivatives at end states
    log_xdot_1st = param.dynfunc_1st_burn(0,log_state(:,1:log_N1st+1),log_control(:,1:log_N1st+1));
    log_xdot_2nd = param.dynfunc_2nd_burn(0,log_state(:,log_N1st+2:end),log_control(:,log_N1st+2:end));
    % states
    C = @(h,x) [1 0 0 0; 0 1 0 0; -3/h^2 -2/h 3/h^2 -1/h; 2/h^3 1/h^2 -2/h^3 1/h^2]*x;
    xspline = @(t,C,h) C(1,floor(t/h)+1)+C(2,floor(t/h)+1).*mod(t,h)+...
        C(3,floor(t/h)+1).*mod(t,h).^2+C(4,floor(t/h)+1).*mod(t,h).^3;
    for i = 1:param.nstate
        % 1st stage phase
        xk = log_state(i,1:log_N1st); xkp1 = log_state(i,2:log_N1st+1);
        xdotk = log_xdot_1st(i,1:end-1); xdotkp1 = log_xdot_1st(i,2:end);
        C_1st = C(log_h1st,[xk; xdotk; xkp1; xdotkp1]);
        init_x((i-1)*N+1:(i-1)*N+N_1st) = xspline(h_1st*(0:N_1st-1),C_1st,log_h1st)';
        % 2nd stage phase
        xk = log_state(i,log_N1st+2:log_N); xkp1 = log_state(i,log_N1st+3:log_N+1);
        xdotk = log_xdot_2nd(i,1:end-1); xdotkp1 = log_xdot_2nd(i,2:end);
        C_2nd = C(log_h2nd,[xk; xdotk; xkp1; xdotkp1]);
        init_x((i-1)*N+N_1st+1:i*N) = xspline(h_2nd*(0:N_2nd-1),C_2nd,log_h2nd)';
    end
    % controls
    ulinear = @(t,uk,ukp1,h) mod(t,h).*(ukp1(floor(t/h)+1)-uk(floor(t/h)+1))./h + uk(floor(t/h)+1);
    for i = 1:param.nctrl
        % 1st stage phase
        uk = log_control(i,1:log_N1st); ukp1 = log_control(i,2:log_N1st+1);
        init_x((param.nstate+i-1)*N+1:(param.nstate+i-1)*N+N_1st) = ulinear(h_1st*(0:N_1st-1),uk,ukp1,log_h1st)';
        % 2nd stage phase
        uk = log_control(i,log_N1st+2:log_N); ukp1 = log_control(i,log_N1st+3:log_N+1);
        init_x((param.nstate+i-1)*N+N_1st+1:(param.nstate+i)*N) = ulinear(h_2nd*(0:N_2nd-1),uk,ukp1,log_h2nd)';
    end
    init_x(ntot*N+1) = h_1st;
    init_x(ntot*N+2) = h_2nd;
    init_x(ntot*N+3) = log_x(log_ntot*log_N+3)/time_scaling;
end
param.axes = plotter(init_x,N_1st,N_2nd,param);

% cost function
cost_func = @(xopt) q_cost(xopt,N_1st,N_2nd,weights,param);
disp("Initial Guess Cost Function: " + num2str(cost_func(init_x)));

% state and control xopt bounds
low_bound = [zeros([N,1]); -20/180*pi*ones([N,1]); zeros([N,1]); param.ms2*ones([N,1]);
    tvc_bounds(1)*ones([N,1]); throtl_bounds(1)*ones([N_1st,1]); throtl_bounds(3)*ones([N_2nd,1]);
    h1st_bounds(1); h2nd_bounds(1); stg_sep_t_bounds(1)];
upp_bound = [sqrt(param.mu/param.earthR)*ones([N,1]); pi/2*ones([N,1]); 3e5*ones([N,1]); param.m0*ones([N,1]);
    tvc_bounds(2)*ones([N,1]); throtl_bounds(2)*ones([N_1st,1]); throtl_bounds(4)*ones([N_2nd,1]);
    h1st_bounds(2); h2nd_bounds(2); stg_sep_t_bounds(2)];

% constraint function
nonl_con_func = @(x) nonl_con_sep(x,N_1st,N_2nd,x0,param);

% optimization solver setup
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
    "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',1e8,'MaxIterations',1e7,...
    EnableFeasibilityMode=true);
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
%     'MaxFunctionEvaluations',1e8,'MaxIterations',1e7,'ConstraintTolerance',1e-6);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
%     "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',1e8,'MaxIterations',1e7,...
%     'ConstraintTolerance',1e-4);
x = fmincon(cost_func,init_x,[],[],[],[],low_bound,upp_bound,nonl_con_func,options);
[csol,ceqsol] = nonl_con_sep(x,N_1st,N_2nd,x0,param);
plotter(x,N_1st,N_2nd,param);

%% FUNCTION - State and Control Plotter
function axes = plotter(x,N_1st,N_2nd,param)
    N = N_1st + N_2nd;
    ntot = param.nstate + param.nctrl;
    t_1st = 0:x(ntot*N+1):x(ntot*N+1)*(N_1st-1);
    t_2nd = x(ntot*N+1)*N_1st:x(ntot*N+2):x(ntot*N+1)*N_1st+x(ntot*N+2)*(N_2nd-1);
    if isfield(param,"axes") axes = param.axes; end
    axes(1) = subplot(2,3,1); plot([t_1st t_2nd]*param.scales.time, x(param.init_ind_ptr(1)+(1:N))*param.scales.speed,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
    axes(2) = subplot(2,3,2); plot([t_1st t_2nd]*param.scales.time, x(param.init_ind_ptr(2)+(1:N))*180/pi,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("Pitch angle (deg)");
    axes(3) = subplot(2,3,3); plot([t_1st t_2nd]*param.scales.time, x(param.init_ind_ptr(3)+(1:N))/1000*param.scales.length,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("Altitude (km)");
    axes(4) = subplot(2,3,4); plot([t_1st t_2nd]*param.scales.time, x(param.init_ind_ptr(4)+(1:N))*param.scales.mass,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");
    axes(5) = subplot(2,3,5); plot([t_1st t_2nd]*param.scales.time, x(param.init_ind_ptr(5)+(1:N))*180/pi,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("Gimbal Angle (deg)");
    axes(6) = subplot(2,3,6); plot([t_1st t_2nd]*param.scales.time, x(param.init_ind_ptr(6)+(1:N))*100,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("Throttle (\%)");
end

%% FUNCTION - Cost function
function [f,g] = q_cost(xopt,N_1st,N_2nd,weights,param)
    N = N_1st + N_2nd;
    ntot = param.nstate + param.nctrl;
    f = -weights.mp_weight*(xopt(param.nstate*N)/param.mp2)^2 + ... % minimize fuel used
        weights.dutvc_weight*sum(diff([xopt(4*N+1:4*N+N_1st)*xopt(ntot*N+1);xopt(4*N+1+N_1st:5*N)*xopt(ntot*N+2)]).^2) + ...
        weights.duT_weight*sum((xopt(5*N+1:6*N)-mean(xopt(5*N+1:6*N))).^2) + ...
        weights.altf_weight*sum(([xopt(2*N+1:2*N+N_1st)*xopt(ntot*N+1);xopt(2*N+1+N_1st:3*N)*xopt(ntot*N+2)]-param.alt_final_ref).^2) + ...
        weights.vf_weight*sum(([xopt(1:N_1st)*xopt(ntot*N+1);xopt(1+N_1st:N)*xopt(ntot*N+2)]-param.v_final_ref).^2);
end

%% FUNCTION - Constraint function
function [c,ceq] = nonl_con_sep(x,N_1st,N_2nd,x0,param)
    persistent func_count;
    if isempty(func_count) func_count = 0; end

    N = N_1st + N_2nd;
    ntot = param.nstate + param.nctrl;
    h_1st = x(ntot*N+1); h_2nd = x(ntot*N+2);
    coast_t = x(ntot*N+3);

    state = [x0 reshape(x(1:param.nstate*N),[N,param.nstate])'];
    control = [[0;1] reshape(x(param.nstate*N+1:ntot*N),[N,param.nctrl])'];
    ceq = nan([param.nstate*N+4,1]);
    init_ind_ptr = param.init_ind_ptr;

    dynfunc_1st_burn = param.dynfunc_1st_burn;
    dynfunc_2nd_burn = param.dynfunc_2nd_burn;

    % Hermite Simpson collocation
    % compute continuous time dynamics at once for end states
    state_1st = state(:,1:N_1st+1); state_2nd = state(:,N_1st+2:end);
    control_1st = control(:,1:N_1st+1); control_2nd = control(:,N_1st+2:end);
    dxdt_1st = dynfunc_1st_burn(0,state_1st,control_1st);
    dxdt_2nd = dynfunc_2nd_burn (0,state_2nd,control_2nd);

    % compute collocation points
    xcol_1st = 1/2*(state_1st(:,1:end-1)+state_1st(:,2:end))-h_1st/8*diff(dxdt_1st,1,2);
    xcol_2nd = 1/2*(state_2nd(:,1:end-1)+state_2nd(:,2:end))-h_2nd/8*diff(dxdt_2nd,1,2);
    xcoldot_1st = 3/2/h_1st*diff(state_1st,1,2)-1/4*(dxdt_1st(:,1:end-1)+dxdt_1st(:,2:end));
    xcoldot_2nd = 3/2/h_2nd*diff(state_2nd,1,2)-1/4*(dxdt_2nd(:,1:end-1)+dxdt_2nd(:,2:end));
    ucol_1st = 1/2*(control_1st(:,1:end-1)+control_1st(:,2:end));
    ucol_2nd = 1/2*(control_2nd(:,1:end-1)+control_2nd(:,2:end));

    % compute continous time dynamics at once for colloction states
    dxdt_col_1st = dynfunc_1st_burn(0,xcol_1st,ucol_1st);
    dxdt_col_2nd = dynfunc_2nd_burn(0,xcol_2nd,ucol_2nd);
    ceq(1:param.nstate*(N-1)) = [reshape((xcoldot_1st-dxdt_col_1st)',[param.nstate*N_1st,1]);
        reshape((xcoldot_2nd-dxdt_col_2nd)',[param.nstate*(N_2nd-1),1])];

    % stage separation constraints from x_N1st to x_N1st+1 reset vehicle mass
    coast_dt = 0.1/param.scales.time;
    dynfunc_stg_sep = @(t,x,u) rocket(t,x,u,param,2);
    stg_sep_record = rk4sim(dynfunc_stg_sep,state(:,N_1st+1)',coast_dt,coast_t,[0 0]);
    ceq(param.nstate*(N-1)+1:param.nstate*N) = ...
        state(:,N_1st+2) - [stg_sep_record.x(end,1:param.nstate-1)';param.m02+param.mplf];
    % ensure initial zero TVC actuation for 2nd stage
    ceq(param.nstate*N+1) = x(init_ind_ptr(param.nstate+1)+N_1st+1,:);

    % final condition constraint
    ceq(param.nstate*N+2) = x(N)-param.v_final_ref; % final error in orbital velocity
    ceq(param.nstate*N+3) = x(2*N)-param.fpa_final_ref; % final error in flight path angle
    ceq(param.nstate*N+4) = x(3*N)-param.alt_final_ref; % final error in orbital altitude

    % propellant usage constraint
    c(1) = param.m0-state(4,N_1st+1)-param.mp1*0.99;
    c(2) = param.m02-state(4,N+1)-param.mp2*0.98;

    % visualization
    if mod(func_count,500) == 0
        param.axes = plotter(x,N_1st,N_2nd,param); drawnow
    end
    func_count = func_count + 1;
end

%% FUNCTION - Generalized RK4 solver for a dynamic equation
function record = rk4sim(dynFunc,init_cond,dt,tf,u_hist)
    t0 = 0; x = init_cond; forces = [0 0];
    t_hist = t0:dt:tf; x_hist = init_cond; force_hist = zeros([1,2]);
    for i = 1:size(t_hist,2)
        t = t_hist(i); x_hist(i,:) = x; force_hist(i,:) = forces;
        if i <= size(u_hist,1) u = u_hist(i,:); else u = [0 0]; end
        [x,forces] = singleRK4(dynFunc,t,dt,x,u);
    end
    record.t = t_hist; record.x = x_hist; record.force = force_hist;
end

%% FUNCTION - single RK4 step
function [next_state,forces] = singleRK4(dynFunc,current_time,dt,current_state,current_control)
    [k1,F1] = dynFunc(current_time,current_state,current_control);
    [k2,F2] = dynFunc(current_time+dt/2,current_state+k1*dt/2,current_control);
    [k3,F3] = dynFunc(current_time+dt/2,current_state+k2*dt/2,current_control);
    [k4,F4] = dynFunc(current_time+dt,current_state+k3*dt,current_control);
    next_state = current_state + 1/6*dt*(k1+2*k2+2*k3+k4);
    forces = 1/6*(F1+2*F2+2*F3+F4);
end

%% FUNCTION - Dynamics of a launch vehicle
function [dxdt,forces] = rocket(t,x,u,param,stage)
% % x represents states as follows
% 1 - v, velocity
% 2 - gamma, flight path angle
% 3 - h, altitude
% 4 - m, total mass

% % u represents control as follows
% 1 - delta TVC angle
% 2 - throttle

% % % single value dynamics
if size(x,1) == 1 || size(x,2) == 1
    v = x(1); fpa = x(2); h = x(3); m = x(4);
    tvc = u(1); throtl = u(2);
    Rp = param.earthR;
    % altitude governed atmosphere
    density = param.rho(h); p = param.P(h); gh = param.g(h);
    % aerodynamic effects
    q = 1/2*density*v^2;
    mach = v/sqrt(1.4*p/density);
    Cd = CD(mach,0); D = q*param.S*Cd;

    if stage == 1
        T = throtl*param.maxT_1st;
        mdot = T/param.Isp1/param.g(0);
    elseif stage == 2
        T = throtl*param.maxT_2nd;
        mdot = T/param.Isp2/param.g(0);
    else
        T = 0; mdot = 0;
    end

    % differential of states
    dxdt = zeros(size(x));
    dxdt(1) = T*cos(tvc)/m-D/m-gh*sin(fpa);
    dxdt(2) = T*sin(tvc)/(m*v)-(gh/v-v/(Rp+h))*cos(fpa);
    dxdt(3) = v*sin(fpa);
    dxdt(4) = -mdot;
    forces = [T D];

else
    v = x(1,:); fpa = x(2,:); h = x(3,:); m = x(4,:);
    tvc = u(1,:); throtl = u(2,:);
    Rp = param.earthR;
    % altitude governed atmosphere
    density = param.rho(h); p = param.P(h); gh = param.g(h);
    % aerodynamic effects
    q = 1/2.*density.*v.^2;
    mach = v./sqrt(1.4.*p./density);
    for i = 1:length(mach) Cd(i) = CD(mach(i),0); end
    D = q.*param.S.*Cd;

    if stage == 1
        T = throtl*param.maxT_1st;
        mdot = T./param.Isp1./param.g(0);
    elseif stage == 2
        T = throtl*param.maxT_2nd;
        mdot = T./param.Isp2./param.g(0);
    else
        T = 0; mdot = 0;
    end

    % differential of states
    dxdt = zeros(size(x));
    dxdt(1,:) = T.*cos(tvc)./m-D./m-gh.*sin(fpa);
    dxdt(2,:) = T.*sin(tvc)./(m.*v)-(gh./v-v./(Rp+h)).*cos(fpa);
    dxdt(3,:) = v.*sin(fpa);
    dxdt(4,:) = -mdot;
    forces = [T D];
end
end

