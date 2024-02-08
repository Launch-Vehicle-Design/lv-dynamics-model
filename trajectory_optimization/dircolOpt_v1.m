clear; clc; close all
addpath("..\parameter_functions")

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

m0 = 1959.8158;
scales.length       = 6357e3;
scales.speed        = sqrt(3.986004418e14/6357e3);
scales.time         = scales.length/scales.speed;
scales.acceleration = scales.speed/scales.time;
scales.mass         = m0;
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
param.h0 = 7640/param.scales.length;
param.g = @(h) 9.80665./(1+h).^2/param.scales.acceleration;
param.rho = @(h) 1.225.*exp(-h/param.h0)/param.scales.density;
param.P = @(h) 101325.*exp(-h/param.h0)/param.scales.pressure;
param.gamma = 1.4;

param.S = pi*(24*0.0254)^2/4/param.scales.area;
param.mPL = 40/param.scales.mass;
% 1st stage
param.m0 = m0/param.scales.mass;
param.ms1 = 104.9236/param.scales.mass;
param.mp1 = 1509.2855/param.scales.mass;
param.Isp1 = 267.6/param.scales.time;
% 2nd stage
param.m02 = 345.6067/param.scales.mass;
param.ms2 = 32.9567/param.scales.mass;
param.mp2 = 266.6499/param.scales.mass;
param.Isp2 = 369.5/param.scales.time;
% payload fairing
param.mplf = 4/param.scales.mass;

param.TtoW_1st = 1.6; param.TtoW_2nd = 0.8;
param.maxT_1st = param.TtoW_1st*param.m0*param.g(0);
param.maxT_2nd = param.TtoW_2nd*param.m02*param.g(0);

init_r = (release_alti+param.earthR)*[cos(release_lat); 0; sin(release_lat)];
init_v = release_velo*[0; -1; 0];
init_cond = [init_r; init_v; param.m0];
% PHASE 1 - uncontrolled drop phase
drop_dt = 0.01/param.scales.time; drop_time = 5/param.scales.time;
drop_dynfun = @(t,x,u) rocket(t,x,u,param,0);
drop_record = rk4sim(drop_dynfun,init_cond,drop_dt,drop_time,zeros([4,1]));

%% Trajectory Optimization - direct trajectory optimization
param.v_final_ref = 7754.1/param.scales.speed;
param.alt_final_ref = 251460/param.scales.length;
param.fpa_final_ref = 0;
param.orb_inc_final_ref = 100/180*pi;
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
weights.vf_weight = 1;
weights.altf_weight = 1;
weights.dutvc_weight = 100;
weights.duT_weight = 10;

% design optimization variable bounds
throtl_bounds = [0.8 1 0.6 1];
h1st_bounds = [0.01 5]/param.scales.time;
h2nd_bounds = [0.01 5]/param.scales.time;
stg_sep_t_bounds = [0.01 60]/param.scales.time;

% design indirect optimization variable constrains
param.tvc_limit = 15/180*pi;

N_1st = 100; N_2nd = 200;
param.N_1st = N_1st; param.N_2nd = N_2nd;
N = N_1st+N_2nd; x0 = drop_record.x(:,end);
% overwrite initial launch attitude
rn0 = norm(x0(1:3)); vn0 = norm(x0(4:6)); v = x0(4:6);
e_rx = x0(1:3)/rn0; e_bz = cross(e_rx,v)/norm(cross(e_rx,v));
pua = -100*pi/180; % pitch up angle
x0(4:6) = cos(pua)*v + sin(pua)*cross(e_bz,v) + (1-cos(pua))*dot(e_bz,v)*e_bz;

% initial guess propagation
param.dynfunc_1st_burn = @(t,x,u) rocket(t,x,u,param,1);
param.dynfunc_2nd_burn = @(t,x,u) rocket(t,x,u,param,2);

%% 1st initial guess - variable range with unified randomness
init_x = zeros([ntot*N+3,1]); 
param.init_ind_ptr = 0:N:ntot*N;
init_x(ntot*N+1) = 1.8/param.scales.time; 
init_x(ntot*N+2) = 2/param.scales.time;
init_x(ntot*N+3) = 0.7/param.scales.time;
init_x(1:param.nstate*N) = reshape([repmat(x0',[N_1st,1]);repmat([x0(1:end-1);param.m02]',[N_2nd,1])],[param.nstate*N,1]);
% init_x(1:param.nstate*N) = reshape([ ...
%     (param.v_final_ref-x0(1))*rand([N,1])+x0(1) ...
%     (param.fpa_final_ref-x0(2))*rand([N,1])+x0(2) ...
%     (param.alt_final_ref-x0(3))*rand([N,1])+x0(3) ...
%     [(param.m0-param.ms1)*rand([N_1st,1])+param.ms1; ...
%     (param.m02-param.mp2)*rand([N_2nd,1])+param.ms2]
%     ],[param.nstate*N,1]);
rand_tvc = rand([3,N]);
init_x(param.nstate*N+1:(param.nstate+3)*N) = reshape((rand_tvc./vecnorm(rand_tvc))',[3*N,1]);
init_x((ntot-1)*N+1:ntot*N) = (throtl_bounds(2)-throtl_bounds(1))*rand([N,1])+throtl_bounds(1);

%% 2nd initial guess - non-optimal educated guess
if user_gues
current_x = x0;
for i = 1:N_1st
    e_v = current_x(4:6)/norm(current_x(4:6));
    e_r = current_x(1:3)/norm(current_x(1:3));
    guess_u = e_r*0.5+2*e_v;
    init_u = [guess_u/norm(guess_u); 0.8];
    dynfunc_1st = param.dynfunc_1st_burn;
    if init_x(ntot*N+1) >= 0.5/param.scales.time
        dt_1st = 0.01/param.scales.time;
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
% % overwrite initial launch attitude
% rn0 = norm(current_x(1:3)); vn0 = norm(current_x(4:6)); v = current_x(4:6);
% e_rx = current_x(1:3)/rn0; e_bz = cross(e_rx,v)/norm(cross(e_rx,v));
% pua = -30*pi/180; % pitch up angle
% current_x(4:6) = cos(pua)*v + sin(pua)*cross(e_bz,v) + (1-cos(pua))*dot(e_bz,v)*e_bz;
current_x(param.nstate) = param.m02;
init_x(param.init_ind_ptr(1:param.nstate)+N_1st+1) = current_x;
init_x(param.init_ind_ptr(param.nstate+1:ntot)+N_1st+1) = init_u;
for i = N_1st+2:N_1st+N_2nd
    e_v = current_x(4:6)/norm(current_x(4:6));
    e_r = current_x(1:3)/norm(current_x(1:3));
    guess_u = e_r*1+1*e_v;
    init_u = [guess_u/norm(guess_u); 1];
    if init_x(ntot*N+2) >= 0.5/param.scales.time
        dt_2nd = 0.01/param.scales.time;
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
file_name = "dircol3d_final1_scaled.mat";
if load_prev && exist(file_name,"file")
    load(file_name,"log_x","log_param");
    state_scaling = ones([7 1]); time_scaling = 1;
    if ~contains(file_name,"scaled")
        state_scaling = [param.scales.length*ones([3 1]); param.scales.speed*ones([3 1]); param.scales.mass];
        time_scaling = param.scales.time;
    end
    init_x = init_guess_interp(x0,log_x,log_param,param,state_scaling,time_scaling);
end

%% 4th initial guess - translate from the 2D solution
file_name = "dircol3d_plfj1st_final1_scaled.mat";
if load_2d && exist(file_name,"file")
    load(file_name,"log_x","log_param");
    state_scaling = [1; 1; 1; 1]; time_scaling = 1;
    if ~contains(file_name,"scaled")
        state_scaling = [param.scales.speed; 1; param.scales.length; param.scales.mass];
        time_scaling = param.scales.time;
    end
    log_N1st = log_param.N_1st; log_N2nd = log_param.N_2nd;
    log_N = log_N1st + log_N2nd; log_ntot = log_param.nstate + log_param.nctrl;
    log_h1st = log_x(log_ntot*log_N+1); log_h2nd = log_x(log_ntot*log_N+2);
    log_ind_ptr = 0:log_N:ntot*log_N;
    % reshape state and control
    state = reshape(log_x(1:log_param.nstate*log_N),[log_N,log_param.nstate])';
    control = reshape(log_x(log_param.nstate*log_N+1:log_ntot*log_N),[log_N,log_param.nctrl])';
    h = state(3,:); v = state(1,:); fpa = state(2,:); m = state(4,:);
    % translate the initial guess optimization variables
    log_init_x = nan([ntot*log_N+3,1]);
    % mass profile does not change
    log_init_x(log_ind_ptr(param.nstate)+(1:log_N)) = m;
    % downrange computation
    ddr = v.*cos(fpa).*[log_h1st*ones([1,log_N1st]) log_h2nd*ones([1,log_N2nd])];
    dr(1) = 0;
    for i = 2:size(ddr,2)
        dr(i) = trapz([ddr(i-1) ddr(i)]) + dr(i-1);
    end
    % position vector
    long = -dr./(param.earthR*cos(release_lat));
    pos = (h+param.earthR).*[
        cos(release_lat)*cos(long);
        cos(release_lat)*sin(long);
        sin(release_lat)*ones(size(h))
        ];
    log_init_x(log_ind_ptr(1)+1:log_ind_ptr(3)+log_N) = reshape(pos',[3*log_N,1]);
    % velocity vector
    v_r = v.*sin(fpa).*[
        cos(release_lat)*cos(long);
        cos(release_lat)*sin(long);
        sin(release_lat)*ones(size(h))
        ];
    v_dr = v.*cos(fpa).*[sin(long); -cos(long); zeros(size(long))];
    log_init_x(log_ind_ptr(4)+1:log_ind_ptr(6)+log_N) = reshape((v_r+v_dr)',[3*log_N,1]);
    % tvc vector
    e_r = [cos(release_lat)*cos(long);cos(release_lat)*sin(long);sin(release_lat)*ones(size(h))];
    e_dr = [sin(long);-cos(long);zeros(size(long))];
    e_bn = cross(e_dr,e_r);
    e_v = (v_r+v_dr)./vecnorm(v_r+v_dr);
    tvc = cos(control(1,:)).*e_v + sin(control(1,:)).*cross(e_bn,e_v) + (1-cos(control(1,:))).*dot(e_bn,e_v).*e_bn;
    log_init_x(log_ind_ptr(param.nstate+1)+1:log_ind_ptr(param.nstate+3)+log_N) = reshape(tvc',[3*log_N,1]);
    % throttle
    log_init_x(log_ind_ptr(param.nstate+param.nctrl)+(1:log_N)) = control(2,:);
    % time translate
    log_init_x(ntot*log_N+(1:3)) = log_x(log_ntot*log_N+(1:3));
    init_x = log_init_x;
end

plotter3d(init_x,N_1st,N_2nd,param);
figure; param.axes = plotter(init_x,N_1st,N_2nd,param);

% cost function
cost_func = @(xopt) q_cost(xopt,N_1st,N_2nd,weights,param);
disp("Initial Guess Cost Function: " + num2str(cost_func(init_x)));

% state and control xopt bounds
low_bound = [-(param.earthR+2*param.alt_final_ref)*ones([3*N,1]); -sqrt(param.mu/param.earthR)*ones([3*N,1]); (param.m0-param.mp1)*ones([N_1st,1]); (param.m02-param.mp2)*ones([N_2nd,1]);
    -1*ones([3*N,1]); throtl_bounds(1)*ones([N_1st,1]); throtl_bounds(3)*ones([N_2nd,1]);
    h1st_bounds(1); h2nd_bounds(1); stg_sep_t_bounds(1)];
upp_bound = [(param.earthR+2*param.alt_final_ref)*ones([3*N,1]); sqrt(param.mu/param.earthR)*ones([3*N,1]); param.m0*ones([N_1st,1]); (param.m02+param.mplf)*ones([N_2nd,1]);
    ones([3*N,1]); throtl_bounds(2)*ones([N_1st,1]); throtl_bounds(4)*ones([N_2nd,1]);
    h1st_bounds(2); h2nd_bounds(2); stg_sep_t_bounds(2)];
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
plotter(x,N_1st,N_2nd,param);

%% FUNCTION - Cost function
function [f,g] = q_cost(xopt,N_1st,N_2nd,weights,param)
    N = N_1st + N_2nd;
    ntot = param.nstate + param.nctrl;
    % state = reshape(xopt(1:param.nstate*N),[N,param.nstate])';
    control = reshape(xopt(param.nstate*N+1:ntot*N),[N,param.nctrl])';
    f = -weights.mp_weight*xopt(param.nstate*N)/param.mp2 + ... % minimize fuel used
        weights.duT_weight*sum((control(end,:)-mean(control(end,:))).^2); % + ...
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
    init_ind_ptr = param.init_ind_ptr;

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
    aug_mass_stg_sep = param.m02+param.mplf;

    % stage separation constraints from x_N1st to x_N1st+1 reset vehicle mass
    coast_dt = 0.1/param.scales.time;
    stg_sep_record = rk4sim(dynfunc_2nd_burn,state(:,N_1st+1),coast_dt,coast_t,[0;0;0;0]);
    ceq(param.nstate*(N-1)+1:param.nstate*N) = state(:,N_1st+2)-[stg_sep_record.x(1:param.nstate-1,end);aug_mass_stg_sep];

    % final condition constraint
    e_H = cross(state(1:3,N),state(4:6,N))/norm(cross(state(1:3,N),state(4:6,N)));
    ceq(param.nstate*N+1) = norm(state(4:6,N))-param.v_final_ref; % final error in orbital velocity
    ceq(param.nstate*N+2) = dot(state(1:3,N),state(4:6,N)); % zero flight path angle
    ceq(param.nstate*N+3) = norm(state(1:3,N))-param.earthR-param.alt_final_ref; % final error in orbital altitude
    ceq(param.nstate*N+4) = acos(dot(e_H,[0;0;1]))-param.orb_inc_final_ref; % orbit plane norm / inclination

    % gimbal sphere angle constraint
    ceq(param.nstate*N+5:(param.nstate+1)*N+4) = vecnorm(control(1:3,2:N+1))-1;

    % propellant usage constraint
    ind_c = 1;
    c(:,ind_c) = param.m0-state(7,N_1st+1)-param.mp1; ind_c = ind_c+1;
    c(:,ind_c) = param.m02-state(7,N+1)-param.mp2; ind_c = ind_c+1;
    % above the ground constraint
    c(:,ind_c:ind_c+N-1) = param.earthR-vecnorm(state(1:3,2:N+1));

    % visualization
    if mod(func_count,500) == 0
        param.axes = plotter(x,N_1st,N_2nd,param); drawnow
    end
    func_count = func_count + 1;
end

%% FUNCTION - Generalized RK4 solver for a dynamic equation
function record = rk4sim(dynFunc,init_cond,dt,tf,u_hist)
    t0 = 0; x = init_cond; forces = zeros([9,1]);
    t_hist = t0:dt:tf; x_hist = init_cond; force_hist = forces;
    for i = 1:size(t_hist,2)
        t = t_hist(i); x_hist(:,i) = x; force_hist(:,i) = forces;
        if i <= size(u_hist,2) u = u_hist(:,i); else u = zeros([4,1]); end
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
function [dxdt,forces] = rocket(t,x,u,param,phase)
% % x represents states as follows 7x1
% 1:3 - r, position of the launch vehicle, iner
% 4:6 - v, velocity of the launch vehicle, iner
% 7 - m, mass of the launch vehicle

% % u represents control as follows 4x1
% 1:3 - TVC direction, iner
% 4 - throttle

% % % single value dynamics
if size(x,1) == 1 || size(x,2) == 1
    r = x(1:3); rn = norm(r);
    v = x(4:6); m = x(7); 
    
    tvc = u(1:3);
    throtl = u(4);

    % derived parameters
    h = rn-param.earthR;
    rh0 = param.rho(h);
    a = sqrt(param.gamma*param.P(h)/rh0);

    % aerodynamic drag
    atmo_v = cross(param.OMEGA,r);
    vrel = v - atmo_v; vreln = norm(vrel);
    e_d = -vrel/vreln;
    Cd = CD(vreln/a,0);
    
    % external forces
    ag = -param.mu/rn^3*r;
    Fd = 1/2*rh0*vreln^2*Cd*param.S*e_d;
    Ft = zeros([3,1]); mdot = 0;
    if phase == 1
        Ft = throtl*param.maxT_1st*tvc;
        mdot = norm(Ft)./param.Isp1./param.g(0);
    elseif phase == 2
        Ft = throtl*param.maxT_2nd*tvc;
        mdot = norm(Ft)./param.Isp2./param.g(0);
    end
    % total acceleration
    atot = ag + (Fd + Ft)/m;
    
    % 1st order state differentials
    dxdt = zeros(size(x));
    dxdt(1:3) = v;
    dxdt(4:6) = atot;
    dxdt(7) = -mdot;
else
    r = x(1:3,:); v = x(4:6,:); m = x(7,:);
    tvc = u(1:3,:); throtl = u(4,:);
    
    % norm array
    rn = vecnorm(r);

    % derived parameters
    h = rn-param.earthR;
    rh0 = param.rho(h);
    a = sqrt(param.gamma.*param.P(h)./rh0);

    % aerodynamic drag
    atmo_v = param.skewOMEGA*r;
    vrel = v - atmo_v; vreln = vecnorm(vrel);
    e_d = -vrel./repmat(vreln,[3,1]);
    mach = vreln./a;
    for i = 1:length(mach) Cd(i) = CD(mach(i),0); end

    % external forces
    ag = -param.mu./rn.^3.*r;
    Fd = 1/2.*rh0.*vreln.^2.*Cd*param.S.*e_d;
    Ft = zeros(size(ag)); mdot = zeros(size(rn));
    if phase == 1
        Ft = throtl.*param.maxT_1st.*tvc;
        mdot = vecnorm(Ft)./param.Isp1./param.g(0);
    elseif phase == 2
        Ft = throtl.*param.maxT_2nd.*tvc;
        mdot = vecnorm(Ft)./param.Isp2./param.g(0);
    end
    % total force
    atot = ag + (Fd + Ft)./m;

    % differential of states
    dxdt = zeros(size(x));
    dxdt(1:3,:) = v;
    dxdt(4:6,:) = atot;
    dxdt(7,:) = -mdot;

end
forces = [ag; Fd; Ft];
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
