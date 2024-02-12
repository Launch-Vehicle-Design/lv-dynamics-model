clear; clc; close all

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

%% Simplified Trajectory Optimization
% % Assumptions
% zero lift generated
% point mass launch vehicle

user_gues = true;
load_prev = true;
load_2d = false;

do_orb_ins_c = false;

mass_mat = [248.214252429984; 27.2989240780479; 174.915328351937;
    1802.36504864883; 93.2490477731307; 1460.90174844571];
if exist("thumper.mat","file")
    load("thumper.mat","optimal","dvdisb");
    mass_mat = [optimal(1:3); optimal(6:8)];
end

scales.length       = 6357e3;
scales.speed        = sqrt(3.986004418e14/6357e3);
scales.time         = scales.length/scales.speed;
scales.acceleration = scales.speed/scales.time;
scales.mass         = mass_mat(4);
scales.force        = scales.mass*scales.acceleration;
scales.area         = scales.length^2;
scales.volume       = scales.area.*scales.length;
scales.density      = scales.mass/scales.volume;
scales.pressure     = scales.force/scales.area;
scales.gravparam    = scales.acceleration*scales.length^2;
param.scales = scales;

release_lat = 80/180*pi;
release_long = -180/180*pi;
release_alti = 40000*0.3048/param.scales.length;
release_velo = 250.786/param.scales.speed;

% general parameter
param.nstate = 7;
param.nctrl = 4;

% environmental parameter
param.mu = 3.986004418e14/param.scales.gravparam;
param.earthR = 6357e3/param.scales.length;
param.OMEGA = [0;0;2*pi/(24*3600)]*param.scales.time;
param.skewOMEGA = [0 -param.OMEGA(3) 0; param.OMEGA(3) 0 0; 0 0 0];
param.gamma = 1.4;
param.g0 = 9.80665/param.scales.acceleration;

param.S = pi*(24*0.0254)^2/4/param.scales.area;
param.mPL = (mass_mat(1)-mass_mat(2)-mass_mat(3))/param.scales.mass;
% 1st stage
param.m0 = mass_mat(4)/param.scales.mass;
param.ms1 = mass_mat(5)/param.scales.mass;
param.mp1 = mass_mat(6)/param.scales.mass;
param.Isp1 = 293.49/param.scales.time;
% 2nd stage
param.m02 = mass_mat(1)/param.scales.mass;
param.ms2 = mass_mat(2)/param.scales.mass;
param.mp2 = mass_mat(3)/param.scales.mass;
param.Isp2 = 369.5/param.scales.time;
% payload fairing
param.mplf = 4/param.scales.mass;

param.TtoW_1st = 1.6; param.TtoW_2nd = 0.8;
param.maxT_1st = param.TtoW_1st*param.m0*param.g0;
param.maxT_2nd = param.TtoW_2nd*param.m02*param.g0;

init_r = (release_alti+param.earthR)*[cos(release_lat); 0; sin(release_lat)];
init_v = release_velo*[0; -1; 0];
init_cond = [init_r; init_v; param.m0];
% PHASE 1 - uncontrolled drop phase
drop_dt = 0.1/param.scales.time; drop_time = 5/param.scales.time;
drop_dynfun = @(t,x,u) rocket(t,x,u,param,0);
drop_record = rk4sim(drop_dynfun,init_cond,drop_dt,drop_time,zeros([4,1]));
param.init_cond = init_cond; param.drop_time = drop_time;

%% Trajectory Optimization - direct trajectory optimization
% final orbit reference
param.v_final_ref = 7754.1/param.scales.speed;
param.alt_final_ref = 251460/param.scales.length;
param.fpa_final_ref = 0;
param.orb_inc_final_ref = 100/180*pi;
param.plf_dropalti_ref = 140e3/param.scales.length;
% constraint reference
param.max_gload = 5;

% xopt represents state and control trajectory of the system
% with total length of 11*N+3(4)
% 1:7N - state from x11...xN1:x51...x5N
% 7N+1:11N - control from u10...u1N-1,u20...u2N-1
% 11N+1:11N+2 - the step length for both 1st and 2nd stage
% 11N+3 - the coast time at stage separation
% 11N+4 - the orbit insertion coast time (optional)

% index helpers
ntot = param.nstate+param.nctrl;

% cost function weighting
weights.mp_weight = 1e2;
weights.fpaf_weight = 10;
weights.vf_weight = 1;
weights.altf_weight = 1;
weights.dutvc_weight = 100;
weights.duT_weight = 10;

% design optimization variable bounds
throtl_bounds = [0.8 1 0.6 1];
h1st_bounds = [0.01 5]/param.scales.time;
h2nd_bounds = [0.01 5]/param.scales.time;
stg_sep_t_bounds = [0.01 60]/param.scales.time;
orb_ins_t_bounds = [0 100]/param.scales.time;

% design indirect optimization variable constrains
param.tvc_limit = 15/180*pi;

N_1st = 100; N_2nd = 200;
param.N_1st = N_1st; param.N_2nd = N_2nd;
N = N_1st+N_2nd; x0 = drop_record.x(:,end);
% overwrite initial launch attitude
rn0 = norm(x0(1:3)); vn0 = norm(x0(4:6)); v = x0(4:6);
e_rx = x0(1:3)/rn0; e_bz = cross(e_rx,v)/norm(cross(e_rx,v));
pua = -90*pi/180; % pitch up angle
x0(4:6) = cos(pua)*v + sin(pua)*cross(e_bz,v) + (1-cos(pua))*dot(e_bz,v)*e_bz;

% initial guess propagation
param.x0 = x0;
param.dynfunc_1st_burn = @(t,x,u) rocket(t,x,u,param,1);
param.dynfunc_2nd_burn = @(t,x,u) rocket(t,x,u,param,2);

%% 1st initial guess - variable range with unified randomness
init_x = zeros([ntot*N+3,1]); 
param.init_ind_ptr = 0:N:ntot*N;
init_x(ntot*N+1) = 1.8/param.scales.time; 
init_x(ntot*N+2) = 2/param.scales.time;
init_x(ntot*N+3) = 0.7/param.scales.time;
if do_orb_ins_c init_x(ntot*N+4) = 0.7/param.scales.time; end
init_x(1:param.nstate*N) = reshape([repmat(x0',[N_1st,1]);repmat([x0(1:end-1);param.m02]',[N_2nd,1])],[param.nstate*N,1]);
rand_tvc = rand([3,N]);
init_x(param.nstate*N+1:(param.nstate+3)*N) = reshape((rand_tvc./vecnorm(rand_tvc))',[3*N,1]);
init_x((ntot-1)*N+1:ntot*N) = (throtl_bounds(2)-throtl_bounds(1))*rand([N,1])+throtl_bounds(1);

%% 2nd initial guess - non-optimal educated guess
if user_gues
current_x = x0;
dt_throtl = 0.1/param.scales.time;
for i = 1:N_1st
    e_v = current_x(4:6)/norm(current_x(4:6));
    e_r = current_x(1:3)/norm(current_x(1:3));
    guess_u = e_r*1+2*e_v;
    init_u = [guess_u/norm(guess_u); 0.8];
    dynfunc_1st = param.dynfunc_1st_burn;
    if init_x(ntot*N+1) >= dt_throtl
        dt_1st = dt_throtl;
        u = repmat(init_u,[1,ceil(init_x(ntot*N+1)/dt_1st)]);
        record = rk4sim(dynfunc_1st,current_x,dt_1st,init_x(ntot*N+1),u);
        current_x = record.x(:,end);
    else
        u = init_u;
        if param.m0 - current_x(4) >= param.mp1 u(2) = 0; end
        current_x = singleRK4(dynfunc_1st,0,init_x(ntot*N+1),current_x,u);
    end
    init_x(param.init_ind_ptr(1:param.nstate)+i) = current_x;
    init_x(param.init_ind_ptr(param.nstate+1:ntot)+i) = init_u;
end
dynfunc_2nd = param.dynfunc_2nd_burn;
stg_sep_record = rk4sim(dynfunc_2nd,current_x,0.01/param.scales.time,init_x(ntot*N+3),[0;0;0;0]);
current_x = stg_sep_record.x(:,end);
init_u = [current_x(4:6)/norm(current_x(4:6)); 1];
current_x(param.nstate) = param.m02;
init_x(param.init_ind_ptr(1:param.nstate)+N_1st+1) = current_x;
init_x(param.init_ind_ptr(param.nstate+1:ntot)+N_1st+1) = init_u;
for i = N_1st+2:N_1st+N_2nd
    e_v = current_x(4:6)/norm(current_x(4:6));
    e_r = current_x(1:3)/norm(current_x(1:3));
    guess_u = e_r*1+2*e_v;
    init_u = [guess_u/norm(guess_u); 1];
    if init_x(ntot*N+2) >= dt_throtl
        dt_2nd = dt_throtl;
        u = repmat(init_u,[1,ceil(init_x(ntot*N+2)/dt_2nd)]);
        record = rk4sim(dynfunc_2nd,current_x,dt_2nd,init_x(ntot*N+2),u);
        current_x = record.x(:,end);
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
file_name = "thumper_straj_cg3dof.mat";
if load_prev && exist(file_name,"file")
    load(file_name,"log_x","log_param");
    state_scaling = ones([7 1]); time_scaling = 1;
    if ~contains(file_name,"scaled") && ~contains(file_name,"straj")
        state_scaling = [param.scales.length*ones([3 1]); param.scales.speed*ones([3 1]); param.scales.mass];
        time_scaling = param.scales.time;
    end
    init_x = init_guess_interp(x0,log_x,log_param,param,state_scaling,time_scaling);
end

plotter3d(init_x,param);
figure; param.axes = plotter(init_x,param);

% cost function
cost_func = @(xopt) q_cost(xopt,N_1st,N_2nd,weights,param);
disp("Initial Guess Cost Function: " + num2str(cost_func(init_x)));

% state and control xopt lower bounds
low_bound = [-(param.earthR+2*param.alt_final_ref)*ones([3*N,1]);
    -sqrt(param.mu/param.earthR)*ones([3*N,1]);
    (param.m0-param.mp1)*ones([N_1st,1]);
    (param.m02-param.mp2)*ones([N_2nd,1]);
    -1*ones([3*N,1]);
    throtl_bounds(1)*ones([N_1st,1]);
    throtl_bounds(3)*ones([N_2nd,1]);
    h1st_bounds(1);
    h2nd_bounds(1);
    stg_sep_t_bounds(1)];

% state and control xopt upper bounds
upp_bound = [(param.earthR+2*param.alt_final_ref)*ones([3*N,1]);
    sqrt(param.mu/param.earthR)*ones([3*N,1]);
    param.m0*ones([N_1st,1]);
    (param.m02+param.mplf)*ones([N_2nd,1]);
    ones([3*N,1]);
    throtl_bounds(2)*ones([N_1st,1]);
    throtl_bounds(4)*ones([N_2nd,1]);
    h1st_bounds(2);
    h2nd_bounds(2);
    stg_sep_t_bounds(2)];

if do_orb_ins_c
    low_bound = [low_bound; orb_ins_t_bounds(1)];
    upp_bound = [upp_bound; orb_ins_t_bounds(2)];
end

ind_lb_violate = init_x <= low_bound;
ind_ub_violate = init_x >= upp_bound;
init_x(ind_lb_violate) = 0.999*low_bound(ind_lb_violate);
init_x(ind_ub_violate) = 0.999*upp_bound(ind_ub_violate);

% constraint function
nonl_con_func = @(x) nonl_con_sep(x,N_1st,N_2nd,x0,param);
nonl_con_func(init_x);

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
plotter(x,param);

%% FUNCTION - Cost function
function [f,g] = q_cost(xopt,N_1st,N_2nd,weights,param)
    N = N_1st + N_2nd;
    % ntot = param.nstate + param.nctrl;
    % state = reshape(xopt(1:param.nstate*N),[N,param.nstate])';
    % control = reshape(xopt(param.nstate*N+1:ntot*N),[N,param.nctrl])';
    f = -weights.mp_weight*xopt(param.nstate*N)/param.mp2; % + ... % minimize fuel used
        % weights.duT_weight*sum((control(end,:)-mean(control(end,:))).^2); % + ...
        % weights.altf_weight*sum([(vecnor/;'m(state(1:3,1:N_1st))'-param.earthR-param.alt_final_ref)*xopt(ntot*N+1);(vecnorm(state(1:3,N_1st+1:N))'-param.earthR-param.alt_final_ref)*xopt(ntot*N+2)].^2) + ...
        % weights.vf_weight*sum([(vecnorm(state(4:6,1:N_1st))'-param.v_final_ref)*xopt(ntot*N+1);(vecnorm(state(4:6,N_1st+1:N))'-param.v_final_ref)*xopt(ntot*N+2)].^2);
    % weights.dutvc_weight*sum(diff([xopt(1*N+1:4*N+N_1st)*xopt(ntot*N+1);xopt(4*N+1+N_1st:5*N)*xopt(ntot*N+2)]).^2) + ...
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
    % set the initial gimbal control to be zero degree
    control = [[x0(4:6)/norm(x0(4:6));1] reshape(x(param.nstate*N+1:ntot*N),[N,param.nctrl])'];
    ceq = nan([(param.nstate+1)*N+4,1]);

    dynfunc_1st_burn = param.dynfunc_1st_burn;
    dynfunc_2nd_burn = param.dynfunc_2nd_burn;

    % Hermite Simpson collocation
    % compute continuous time dynamics at once for end states
    state_1st = state(:,1:N_1st+1); state_2nd = state(:,N_1st+2:end);
    control_1st = control(:,1:N_1st+1); control_2nd = control(:,N_1st+2:end);
    dxdt_1st = dynfunc_1st_burn(0,state_1st,control_1st);
    dxdt_2nd = dynfunc_2nd_burn (0,state_2nd,control_2nd);

    % locate the payload fairing drop point
    ind_height_abv = find(vecnorm(state(1:3,:))>=param.plf_dropalti_ref);
    ind_jet = ind_height_abv(1); jet_2nd = 2+N_1st-ind_jet<0;

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
    aug_mass_stg_sep = param.m02;
    ceq(1:param.nstate*(N-1)) = [reshape((xcoldot_1st-dxdt_col_1st)',[param.nstate*N_1st,1]);
        reshape((xcoldot_2nd-dxdt_col_2nd)',[param.nstate*(N_2nd-1),1])];
    % if jet_2nd
    %     ind_jet2nd = ind_jet-(1+N_1st);
    %     aug_mass_stg_sep = param.m02+param.mplf;
    %     ceq(1:param.nstate*(N-2)) = [reshape((xcoldot_1st-dxdt_col_1st)',[param.nstate*N_1st,1]);
    %         reshape((xcoldot_2nd(:,1:ind_jet2nd-1)-dxdt_col_2nd(:,1:ind_jet2nd-1))',[param.nstate*(ind_jet2nd-1),1]);
    %         reshape((xcoldot_2nd(:,ind_jet2nd+1:end)-dxdt_col_2nd(:,ind_jet2nd+1:end))',[param.nstate*(N_2nd-ind_jet2nd-1),1])];
    %     ceq(param.nstate*(N-2)+1:param.nstate*(N-1)-1) = xcoldot_2nd(1:end-1,ind_jet2nd)-dxdt_col_2nd(1:end-1,ind_jet2nd);
    %     ceq(param.nstate*(N-1)) = state_2nd(end,ind_jet2nd+1)+param.mplf-dxdt_col_2nd(end,ind_jet2nd)*h_2nd-state_2nd(end,ind_jet2nd);
    % end
    % aug_mass_stg_sep = param.m02+param.mplf;

    % stage separation constraints from x_N1st to x_N1st+1 reset vehicle mass
    coast_dt = min([10/param.scales.time coast_t]);
    stg_sep_record = rk4sim(dynfunc_2nd_burn,state(:,N_1st+1),coast_dt,coast_t,[0;0;0;0]);
    ceq(param.nstate*(N-1)+1:param.nstate*N) = state(:,N_1st+2)-[stg_sep_record.x(1:param.nstate-1,end);aug_mass_stg_sep];

    % 2nd stage MECO coast to orbit insertion
    final_state = state(:,N);
    if size(x,2) == ntot*N+4
        orb_ins_t = x(ntot*N+4);
        orb_ins_dt = min([10/param.scales.time orb_ins_t]);
        orb_ins_record = rk4sim(dynfunc_2nd_burn,state(:,N),orb_ins_dt,orb_ins_t,[0;0;0;0]);
        final_state = orb_ins_record.x(:,end);
    end

    % final condition constraint
    e_H = cross(final_state(1:3),final_state(4:6))/norm(cross(final_state(1:3),final_state(4:6)));
    ceq(param.nstate*N+1) = norm(final_state(4:6))-param.v_final_ref; % final error in orbital velocity
    ceq(param.nstate*N+2) = dot(final_state(1:3),final_state(4:6)); % zero flight path angle
    ceq(param.nstate*N+3) = norm(final_state(1:3))-param.earthR-param.alt_final_ref; % final error in orbital altitude
    ceq(param.nstate*N+4) = acos(dot(e_H,[0;0;1]))-param.orb_inc_final_ref; % orbit plane norm / inclination

    % gimbal sphere angle constraint
    ceq(param.nstate*N+5:(param.nstate+1)*N+4) = vecnorm(control(1:3,2:N+1))-1;

    % propellant usage constraint
    ind_c = 1;
    c(:,ind_c) = param.m0-state(7,N_1st+1)-param.mp1; ind_c = ind_c+1;
    c(:,ind_c) = param.m02-state(7,N+1)-param.mp2; ind_c = ind_c+1;

    % above the ground constraint
    c(:,ind_c:ind_c+N-1) = param.earthR-vecnorm(state(1:3,2:N+1)); ind_c = ind_c+N;
    
    % payload g-load constraint
    c(:,ind_c:ind_c+N-2) = [vecnorm(dxdt_col_1st(4:6,:)) vecnorm(dxdt_col_2nd(4:6,:))]-param.g0*param.max_gload;

    % visualization
    if mod(func_count,500) == 0
        param.axes = plotter(x,param); drawnow
    end
    func_count = func_count + 1;
end

%% FUNCTION - initial guess interpolation from previous result
function init_x = init_guess_interp(x0,log_x,log_param,param,state_scaling,time_scaling)
    N_1st = param.N_1st; N_2nd = param.N_2nd; N = N_1st + N_2nd;
    ntot = param.nstate + param.nctrl;
    log_N1st = log_param.N_1st; log_N2nd = log_param.N_2nd;
    log_N = log_N1st + log_N2nd; log_ntot = log_param.nstate + log_param.nctrl;
    log_h1st = log_x(log_ntot*log_N+1)/time_scaling;
    log_h2nd = log_x(log_ntot*log_N+2)/time_scaling;
    log_state = [x0 reshape(log_x(1:log_param.nstate*log_N),[log_N,log_param.nstate])'./state_scaling];
    log_control = [[x0(4:6)/norm(4:6);1] reshape(log_x(log_param.nstate*log_N+1:log_ntot*log_N),[log_N,log_param.nctrl])'];
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
