clear; clc; close all
addpath("..\parameter_functions")

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

%% Simplified Trajectory Optimization
% % Assumptions
% zero angle of attack
% zero lift generated
% point mass launch vehicle

user_gues = true;
load_prev = false;

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

release_lat = 80/180*pi; 
release_long = -180/180*pi;
release_alti = 40000*0.3048/param.scales.length;
release_velo = 250.786/param.scales.speed;

% general parameter
param.nstate = 7;
param.nctrl = 3;

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
param.m0 = 2241.1532/param.scales.mass;
param.ms1 = 125.2035/param.scales.mass;
param.mp1 = 1801.0036/param.scales.mass;
param.Isp1 = 276/param.scales.time;
% 2nd stage
param.m02 = 314.9462/param.scales.mass;
param.ms2 = 29.1843/param.scales.mass;
param.mp2 = 248.7618/param.scales.mass;
param.Isp2 = 369.5/param.scales.time;

param.TtoW_1st = 1.6; param.TtoW_2nd = 0.8;
param.maxT_1st = param.TtoW_1st*param.m0*param.g(0);
param.maxT_2nd = param.TtoW_2nd*param.m02*param.g(0);

init_r = (release_alti+param.earthR)*[cos(release_lat); 0; sin(release_lat)];
init_v = release_velo*[0; -1; 0];
init_cond = [init_r; init_v; param.m0];
% PHASE 1 - uncontrolled drop phase
drop_dt = 0.01/param.scales.time; drop_time = 5/param.scales.time;
drop_dynfun = @(t,x,u) rocket(t,x,u,param,0);
drop_record = rk4sim(drop_dynfun,init_cond,drop_dt,drop_time,zeros([3,1]));

%% Trajectory Optimization - direct trajectory optimization
param.v_final_ref = 7754.1/param.scales.speed;
param.alt_final_ref = 251460/param.scales.length;
param.fpa_final_ref = 0;
param.orb_inc_final_ref = 100/180*pi;
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
tvc_bounds = [-15 15]/180*pi;

N_1st = 100; N_2nd = 100;
N = N_1st+N_2nd; x0 = drop_record.x(:,end);
% overwrite initial launch attitude
rn0 = norm(x0(1:3)); v = x0(4:6);
e_rx = x0(1:3)/rn0; e_bz = cross(e_rx,v)/norm(cross(e_rx,v));
pua = -90*pi/180; % pitch up angle
x0(4:6) = cos(pua)*v + sin(pua)*cross(e_bz,v) + (1-cos(pua))*dot(e_bz,v)*e_bz;

% initial guess propagation
param.dynfunc_1st_burn = @(t,x,u) rocket(t,x,u,param,1);
param.dynfunc_2nd_burn = @(t,x,u) rocket(t,x,u,param,2);

%% 1st initial guess - variable range with unified randomness
init_x = zeros([ntot*N+3,1]); 
param.init_ind_ptr = 0:N:ntot*N;
init_x(ntot*N+1) = 1.6/param.scales.time; 
init_x(ntot*N+2) = 2/param.scales.time;
init_x(ntot*N+3) = 0.5/param.scales.time;
init_x(1:param.nstate*N) = reshape([repmat(x0',[N_1st,1]);repmat([x0(1:end-1);param.m02]',[N_2nd,1])],[param.nstate*N,1]);
% init_x(1:param.nstate*N) = reshape([ ...
%     (param.v_final_ref-x0(1))*rand([N,1])+x0(1) ...
%     (param.fpa_final_ref-x0(2))*rand([N,1])+x0(2) ...
%     (param.alt_final_ref-x0(3))*rand([N,1])+x0(3) ...
%     [(param.m0-param.ms1)*rand([N_1st,1])+param.ms1; ...
%     (param.m02-param.mp2)*rand([N_2nd,1])+param.ms2]
%     ],[param.nstate*N,1]);
rand_tvc = [(diff(tvc_bounds)+tvc_bounds(1))*rand([N,1]); 2*pi*rand([N,1])];
init_x(param.nstate*N+1:(param.nstate+2)*N) = rand_tvc;
init_x((ntot-1)*N+1:ntot*N) = (throtl_bounds(2)-throtl_bounds(1))*rand([N,1])+throtl_bounds(1);

%% 2nd initial guess - non-optimal educated guess
if user_gues
current_x = x0;
dynfunc_1st = param.dynfunc_1st_burn;
for i = 1:N_1st
    rn = norm(current_x(1:3)); vn = norm(current_x(4:6));
    pitch_ang = pi/2-acos((dot(current_x(1:3),current_x(4:6)))./rn./vn);
    init_u = [max(0,min(15/180*pi,-(pitch_ang-60*pi/180))); sign(60*pi/180-pitch_ang)*pi/2; 0.8];
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
stg_sep_record = rk4sim(dynfunc_2nd,current_x,0.01/param.scales.time,init_x(ntot*N+3),[0;0;0]);
current_x = stg_sep_record.x(:,end);
current_x(param.nstate) = param.m02;
init_u = [0;0;1];
init_x(param.init_ind_ptr(1:param.nstate)+N_1st+1) = current_x;
init_x(param.init_ind_ptr(param.nstate+1:ntot)+N_1st+1) = init_u;
for i = N_1st+2:N_1st+N_2nd
    rn = norm(current_x(1:3)); vn = norm(current_x(4:6));
    pitch_ang = pi/2-acos((dot(current_x(1:3),current_x(4:6)))./rn./vn);
    init_u = [max(0,min(15/180*pi,-(pitch_ang-60*pi/180))); sign(60*pi/180-pitch_ang)*pi/2; 1];
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
file_name = "dircol_sqp_final1_scaled.mat";
if load_prev && exist(file_name,"file")
    load(file_name,"log_x","log_param");
    state_scaling = [1; 1; 1; 1]; time_scaling = 1;
    if ~contains(file_name,"scaled")
        state_scaling = [param.scales.speed; 1; param.scales.length; param.scales.mass];
        time_scaling = param.scales.time;
    end
    % translate 2D to 3D
    state = reshape(x(1:log_param.nstate*N),[N,log_param.nstate])';
    control = reshape(x(log_param.nstate*N+1:ntot*N),[N,log_param.nctrl])';
    v0 = state(1,1); r0 = state(3,1); gam0 = state(2,1);
    init_state = []

    init_x(ntot*N+1) = h_1st;
    init_x(ntot*N+2) = h_2nd;
    init_x(ntot*N+3) = log_x(log_ntot*log_N+3)/time_scaling;
end

plotter3d(init_x,N_1st,N_2nd,param);
figure; param.axes = plotter(init_x,N_1st,N_2nd,param);

% cost function
cost_func = @(xopt) q_cost(xopt,N_1st,N_2nd,weights,param);
disp("Initial Guess Cost Function: " + num2str(cost_func(init_x)));

% state and control xopt bounds
low_bound = [-(param.earthR+2*param.alt_final_ref)*ones([3*N,1]); -sqrt(param.mu/param.earthR)*ones([3*N,1]); (param.m0-param.mp1)*ones([N_1st,1]); (param.m02-param.mp2)*ones([N_1st,1]);
    tvc_bounds(1)*ones([N,1]); zeros([N,1]); throtl_bounds(1)*ones([N_1st,1]); throtl_bounds(3)*ones([N_2nd,1]);
    h1st_bounds(1); h2nd_bounds(1); stg_sep_t_bounds(1)];
upp_bound = [(param.earthR+2*param.alt_final_ref)*ones([3*N,1]); sqrt(param.mu/param.earthR)*ones([3*N,1]); param.m0*ones([N_1st,1]); param.m02*ones([N_2nd,1]);
    tvc_bounds(2)*ones([N,1]); 2*pi*ones([N,1]); throtl_bounds(2)*ones([N_1st,1]); throtl_bounds(4)*ones([N_2nd,1]);
    h1st_bounds(2); h2nd_bounds(2); stg_sep_t_bounds(2)];

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
    state = reshape(xopt(1:param.nstate*N),[N,param.nstate])';
    control = reshape(xopt(param.nstate*N+1:ntot*N),[N,param.nctrl])';
    f = -weights.mp_weight*(xopt(param.nstate*N)/param.mp2)^2 + ... % minimize fuel used
        weights.duT_weight*sum((control(end,:)-mean(control(end,:))).^2) + ...
        weights.altf_weight*sum(([vecnorm(state(1:3,1:N_1st))'*xopt(ntot*N+1);vecnorm(state(1:3,N_1st+1:N))'*xopt(ntot*N+2)]-param.earthR-param.alt_final_ref).^2) + ...
        weights.vf_weight*sum(([vecnorm(state(4:6,1:N_1st))'*xopt(ntot*N+1);vecnorm(state(4:6,N_1st+1:N))'*xopt(ntot*N+2)]-param.v_final_ref).^2);
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
    control = [[0;0;1] reshape(x(param.nstate*N+1:ntot*N),[N,param.nctrl])'];
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
    stage_sep_x = [state(1:param.nstate-1,N_1st+1); param.m02];
    stg_sep_record = rk4sim(dynfunc_2nd_burn,stage_sep_x,coast_dt,coast_t,[0;0;0]);
    ceq(param.nstate*(N-1)+1:param.nstate*N) = state(:,N_1st+2) - stg_sep_record.x(:,end);

    % final condition constraint
    e_H = cross(state(1:3,N),state(4:6,N))/norm(cross(state(1:3,N),state(4:6,N)));
    ceq(param.nstate*N+1) = norm(state(4:6,N))-param.v_final_ref; % final error in orbital velocity
    ceq(param.nstate*N+2) = dot(state(1:3,N),state(4:6,N)); % zero flight path angle
    ceq(param.nstate*N+3) = norm(state(1:3,N))-param.earthR-param.alt_final_ref; % final error in orbital altitude
    ceq(param.nstate*N+4) = acos(dot(e_H,[0;0;1]))-param.orb_inc_final_ref; % orbit plane norm / inclination

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
        if i <= size(u_hist,2) u = u_hist(:,i); else u = zeros([3,1]); end
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
% 1:2 - TVC direction angles
% 3 - throttle

% % % single value dynamics
if size(x,1) == 1 || size(x,2) == 1
    r = x(1:3); rn = norm(r);
    v = x(4:6); vn = norm(v); 
    m = x(7);
    
    tvc = u(1:2);
    throtl = u(3);

    % derived parameters
    h = rn-param.earthR;
    rh0 = param.rho(h);
    a = sqrt(param.gamma*param.P(h)/rh0);

    % aerodynamic drag
    atmo_v = cross(param.OMEGA,r);
    vrel = v - atmo_v; vreln = norm(vrel);
    e_d = -vrel/vreln;
    Cd = CD(vreln/a,0);

    % TVC direction
    e_down = -r/rn;
    east = cross(e_down,[0;0;1]); e_east = east/norm(east);
    e_south = cross(e_down,e_east);
    e_vz = v/vn;
    bx = cross(e_vz,e_down);
    if norm(bx) <= 1e-6
        e_bx = e_south;
    else
        e_bx = bx/norm(bx);
    end
    e_ny = cross(e_vz,e_bx)/norm(cross(e_vz,e_bx));
    e_tvc = sin(tvc(1))*cos(tvc(2))*e_bx + sin(tvc(1))*sin(tvc(2))*e_ny + cos(tvc(1))*e_vz;
    
    % external forces
    ag = param.mu/rn^2*e_down;
    Fd = 1/2*rh0*vreln^2*Cd*param.S*e_d;
    Ft = zeros([3,1]); mdot = 0;
    if phase == 1
        Ft = throtl*param.maxT_1st*e_tvc;
        mdot = norm(Ft)./param.Isp1./param.g(0);
    elseif phase == 2
        Ft = throtl*param.maxT_2nd*e_tvc;
        mdot = norm(Ft)./param.Isp2./param.g(0);
    end
    % total force
    atot = ag + (Fd + Ft)/m;
    
    % 1st order state differentials
    dxdt = zeros(size(x));
    dxdt(1:3) = v;
    dxdt(4:6) = atot;
    dxdt(7) = -mdot;
else
    r = x(1:3,:); v = x(4:6,:); m = x(7,:);
    tvc = u(1:2,:); throtl = u(3,:);
    
    % norm array
    rn = vecnorm(r); vn = vecnorm(v);

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

    % TVC direction
    e_down = -r./rn;
    east = cross(e_down,repmat([0;0;1],[1,length(r)])); e_east = east./vecnorm(east);
    e_south = cross(e_down,e_east);
    e_vz = v./vn;
    bx = cross(e_vz,e_down);
    e_bx = bx./vecnorm(bx);
    verti_ind = vecnorm(bx) <= 1e-6;
    e_bx(:,verti_ind) = e_south(:,verti_ind);
    ny = cross(e_vz,e_bx); e_ny = ny./vecnorm(ny);
    e_tvc = sin(tvc(1,:)).*cos(tvc(2,:)).*e_bx + sin(tvc(1)).*sin(tvc(2)).*e_ny + cos(tvc(1)).*e_vz;

    % external forces
    ag = param.mu./rn.^2.*e_down;
    Fd = 1/2.*rh0.*vreln.^2.*Cd*param.S.*e_d;
    Ft = zeros(size(ag)); mdot = zeros(size(rn));
    if phase == 1
        Ft = throtl.*param.maxT_1st.*e_tvc;
        mdot = vecnorm(Ft)./param.Isp1./param.g(0);
    elseif phase == 2
        Ft = throtl.*param.maxT_2nd.*e_tvc;
        mdot = vecnorm(Ft)./param.Isp2./param.g(0);
    end
    % total acceleration
    atot = ag + (Fd + Ft)./m;

    % differential of states
    dxdt = zeros(size(x));
    dxdt(1:3,:) = v;
    dxdt(4:6,:) = atot;
    dxdt(7,:) = -mdot;

end
forces = [ag; Fd; Ft];
end

