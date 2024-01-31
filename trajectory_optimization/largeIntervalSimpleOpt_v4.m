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

release_alti = 40000*0.3048;
release_velo = 250.786;

% general parameter
param.nstate = 4;
param.nctrl = 2;

% environmental parameter
param.mu = 3.986004418e14;
param.earthR = 6357e3;
param.g = @(h) g(h);
param.rho = @(h) rho(h);
param.P = @(h) P(h);
param.gamma = 1.4;

param.lPLF = 26*0.0254;
param.dia = 24*0.0254;
param.S = pi*param.dia^2/4;
param.mPL = 40;
% 1st stage
param.m0 = 2241.1532;
param.ms1 = 125.2035;
param.mp1 = 1801.0036;
param.Isp1 = 276;
param.l1 = 155.7*0.0254;
% 2nd stage
param.m02 = 314.9462;
param.ms2 = 29.1843;
param.mp2 = 248.7618;
param.Isp2 = 369.5;
param.l2 = 55.6*0.0254;

param.TtoW_1st = 1.6; param.TtoW_2nd = 0.8;
param.maxT_1st = param.TtoW_1st*param.m0*9.80665;
param.maxT_2nd = param.TtoW_2nd*param.m02*9.80665;

param.scale.velo = sqrt(param.mu/param.earthR);
param.scale.alti = 100e3;
param.scale.fpa = pi/180;
param.scale.mass = 1e3;

init_cond = [release_velo 0 release_alti param.m0];
% PHASE 1 - uncontrolled drop phase
drop_dt = 0.01; drop_time = 5;
drop_dynfun = @(t,x,u) rocket(t,x,u,param,0);
drop_record = rk4sim(drop_dynfun,init_cond,drop_dt,drop_time,[0 0]);

%% Trajectory Optimization - direct trajectory optimization
param.v_final_ref = 7754.1;
param.alt_final_ref = 251460;
param.fpa_final_ref = 0;
% xopt represents state and control trajectory of the system
% with total length of 7*N+1
% 1:4N - state from x11...xN1:x51...x5N
% 4N+1:6N - control from u10...u1N-1,u20...u2N-1
% 6N+1:6N+2 - the step length for both 1st and 2nd stage
% 6N+3 - the coast time at stage separation

% index helpers
ntot = param.nstate+param.nctrl;

% cost function weighting
mp_weight = 1;
fpaf_weight = 10;
vf_weight = 10;
altf_weight = 10;

% design optimal variable bounds
tvc_bounds = [-15 15]/180*pi;
throtl_bounds = [0.4 1.01];
h1st_bounds = [0.5 30];
h2nd_bounds = [0.5 50];
stg_sep_t_bounds = [0.1 100];

N_1st = 10; N_2nd = 10;
N = N_1st+N_2nd; x0 = drop_record.x(end,:)';
x0(2) = 80*pi/180;

% initial guess propagation
init_x = zeros([ntot*N+3,1]); 
param.init_ind_ptr = 0:N:ntot*N;
init_x(ntot*N+1) = 18; 
init_x(ntot*N+2) = 32; 
init_x(ntot*N+3) = 5;
% init_x(1:5*N) = reshape(repmat([param.v_final_ref param.fpa_final_ref param.alt_final_ref param.m02-param.mp2],[N,1]),[5*N,1]);
% init_x(1:5*N) = reshape([ ...
%     (param.v_final_ref-param.x0(1))*rand([N,1])+param.x0(1) ...
%     (param.fpa_final_ref-param.x0(2))*rand([N,1])+param.x0(2) ...
%     (param.alt_final_ref-param.x0(3))*rand([N,1])+param.x0(3) ...
%     [(param.m0-param.ms1)*rand([N_1st,1])+param.ms1; ...
%     (param.m02-param.mp2)*rand([N_2nd,1])+param.ms2]
%     ],[5*N,1]);
% init_x(5*N+1:6*N) = (tvc_bounds(2)-tvc_bounds(1))*rand([N,1])+tvc_bounds(1);
% init_x(6*N+1:7*N) = (throtl_bounds(2)-throtl_bounds(1))*rand([N,1])+throtl_bounds(1);

% % % 2nd initial guess
% current_x = x0;
% for i = 1:N_1st
%     init_u = [max(min(-1*(current_x(2)-50/180*pi),tvc_bounds(2)*0.99),tvc_bounds(1)*0.99) 0.8];
%     dynfunc_1st = @(t,x,u) rocket(t,x,u,param,1);
%     if init_x(ntot*N+1) >= 0.5
%         dt_1st = 0.01;
%         u = repmat(init_u,[ceil(init_x(ntot*N+1)/dt_1st),1]);
%         record = rk4sim(dynfunc_1st,current_x',dt_1st,init_x(ntot*N+1),u);
%         current_x = record.x(end,:)';
%     else
%         u = init_u;
%         if param.m0 - current_x(4) >= param.mp1 u(2) = 0; end
%         current_x = singleRK4(dynfunc_1st,0,init_x(ntot*N+1),current_x,u);
%     end
%     init_x(param.init_ind_ptr(1:param.nstate)+i) = current_x;
%     init_x(param.init_ind_ptr(param.nstate+1:ntot)+i) = init_u;
% end
% current_x(param.nstate) = param.m02;
% dynfunc_2nd = @(t,x,u) rocket(t,x,u,param,2);
% stg_sep_record = rk4sim(dynfunc_2nd,current_x',0.1,init_x(ntot*N+3),[0 0]);
% current_x = stg_sep_record.x(end,:)';
% init_u = [0 1];
% init_x(param.init_ind_ptr(1:param.nstate)+N_1st+1) = current_x;
% init_x(param.init_ind_ptr(param.nstate+1:ntot)+N_1st+1) = init_u;
% for i = N_1st+2:N_1st+N_2nd
%     init_u = [max(min(-1*(current_x(2)-30/180*pi),tvc_bounds(2)*0.99),tvc_bounds(1)*0.99) 1];
%     if init_x(ntot*N+2) >= 0.5
%         dt_2nd = 0.01;
%         u = repmat(init_u,[ceil(init_x(ntot*N+2)/dt_2nd),1]);
%         record = rk4sim(dynfunc_2nd,current_x',dt_2nd,init_x(ntot*N+2),u);
%         current_x = record.x(end,:)';
%     else
%         u = init_u;
%         if param.m02 - current_x(4) >= param.mp2 u(2) = 0; end
%         current_x = singleRK4(dynfunc_2nd,0,init_x(ntot*N+2),current_x,u);
%     end
%     init_x(param.init_ind_ptr(1:param.nstate)+i) = current_x;
%     init_x(param.init_ind_ptr(param.nstate+1:ntot)+i) = init_u;
% end

% % % 3rd initial guess
if exist("progress.mat","file")
    load("progress.mat","log_x");
    init_x = log_x;
end

plotter(init_x,drop_time,N_1st,N_2nd,param);

% cost function
cost_func = @(xopt) -mp_weight*xopt(ntot*N)^2 + ... % minimize fuel used
    vf_weight*((xopt(N)-param.v_final_ref)/param.scale.velo)^2 + ... % final error in orbital velocity
    fpaf_weight*((xopt(2*N)-param.fpa_final_ref)/param.scale.fpa)^2 + ... % final error in flight path angle
    altf_weight*((xopt(3*N)-param.alt_final_ref)/param.scale.alti)^2 + ... % final error in orbital altitude
    sum(((xopt(2*N+1:3*N)-param.alt_final_ref)/param.scale.alti).^2) + ...
    sum(((xopt(1:N)-param.v_final_ref)/param.scale.velo).^2);
disp("Initial Guess Cost Function: " + num2str(cost_func(init_x)));

% state and control xopt bounds
low_bound = [-inf*ones([param.nstate*N,1]);
    tvc_bounds(1)*ones([N,1]); throtl_bounds(1)*ones([N,1]);
    h1st_bounds(1); h2nd_bounds(1); stg_sep_t_bounds(1)];
upp_bound = [inf*ones([param.nstate*N,1]);
    tvc_bounds(2)*ones([N,1]); throtl_bounds(2)*ones([N,1]);
    h1st_bounds(2); h2nd_bounds(2); stg_sep_t_bounds(2)];

% constraint function
nonl_con_func = @(x) nonl_con_sep(x,N_1st,N_2nd,x0,param);

% optimization solver setup
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
    "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',1e8,'MaxIterations',1e7,...
    EnableFeasibilityMode=true);
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
    'MaxFunctionEvaluations',1e8,'MaxIterations',1e7);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
%     "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',1e8,'MaxIterations',1e7);
x = fmincon(cost_func,init_x,[],[],[],[],low_bound,upp_bound,nonl_con_func,options);
[csol,ceqsol] = nonl_con_sep(x,N_1st,N_2nd,x0,param);
plotter(x,drop_time,N_1st,N_2nd,param);

%% FUNCTION - State and Control Plotter
function plotter(x,drop_time,N_1st,N_2nd,param)
N = N_1st + N_2nd;
ntot = param.nstate + param.nctrl;
t_1st = drop_time:x(ntot*N+1):drop_time+x(ntot*N+1)*(N_1st-1);
t_2nd = drop_time+x(ntot*N+1)*N_1st:x(ntot*N+2):drop_time+x(ntot*N+1)*N_1st+x(ntot*N+2)*(N_2nd-1);
figure; subplot(2,3,1); plot([t_1st t_2nd], x(param.init_ind_ptr(1)+(1:N)),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
subplot(2,3,2); plot([t_1st t_2nd], x(param.init_ind_ptr(2)+(1:N))*180/pi,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Pitch angle (deg)");
subplot(2,3,3); plot([t_1st t_2nd], x(param.init_ind_ptr(3)+(1:N))/1000,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Altitude (km)");
subplot(2,3,5); plot([t_1st t_2nd], x(param.init_ind_ptr(4)+(1:N)),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");
subplot(2,3,6); plot([t_1st t_2nd], x(param.init_ind_ptr(5)+(1:N))*180/pi,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Gimbal Angle (deg)");
hold on; plot([t_1st t_2nd], x(param.init_ind_ptr(6)+(1:N)),"k","LineWidth",1.2);
end

%% FUNCTION - Constraint function
function [c,ceq] = nonl_con_sep(x,N_1st,N_2nd,x0,param)
    N = N_1st + N_2nd;
    ntot = param.nstate + param.nctrl;
    h_1st = x(ntot*N+1); h_2nd = x(ntot*N+2);
    coast_t = x(ntot*N+3);

    state = [x0 reshape(x(1:param.nstate*N),[N,param.nstate])'];
    ceq = nan([param.nstate*N+4,1]);
    init_ind_ptr = param.init_ind_ptr;

    dynfunc_1st_burn = @(t,x,u) rocket(t,x,u,param,1);
    dynfunc_2nd_burn = @(t,x,u) rocket(t,x,u,param,2);

    % equality constraint from dynamics of the 1st stage burn
    % from x0 to x_N1st and u0 to u_N1st-1
    % total of N1st+1 states and N1st controls
    for i = 1:N_1st
        current_x = state(:,i);
        current_u = x(init_ind_ptr(param.nstate+1:ntot)+i,:);
        next_x = state(:,i+1);
        if h_1st >= 0.5
            dt_1st = 0.1;
            u = repmat(current_u',[ceil(h_1st/dt_1st),1]);
            step_record = rk4sim(dynfunc_1st_burn,current_x',dt_1st,h_1st,u);
            ceq(param.nstate*(i-1)+1:param.nstate*i) = next_x - step_record.x(end,:)';
        else
            ceq(param.nstate*(i-1)+1:param.nstate*i) = next_x - singleRK4(dynfunc_1st_burn,0,h_1st,current_x,current_u);
        end
    end

    % stage separation constraints from x_N1st to x_N1st+1 reset vehicle mass
    coast_dt = 0.1;
    stage_sep_x = [state(1:param.nstate-1,N_1st+1); param.m02];
    dynfunc_stg_sep = @(t,x,u) rocket(t,x,u,param,2);
    stg_sep_record = rk4sim(dynfunc_stg_sep,stage_sep_x',coast_dt,coast_t,[0 0]);
    ceq(param.nstate*N_1st+1:param.nstate*(N_1st+1)) = state(:,N_1st+2) - stg_sep_record.x(end,:)';
    % ensure initial zero TVC actuation for 2nd stage
    ceq(param.nstate*N+1) = x(init_ind_ptr(param.nstate+1)+N_1st+1,:);
    
    % equality constraint from dynamics of the 2nd stage burn
    % from x_N1st+1 to x_N and uN1st+1 to u_N-1
    % total of N2nd states and N2nd-1 controls
    for i = N_1st+2:N_1st+N_2nd
        current_x = state(:,i);
        current_u = x(init_ind_ptr(param.nstate+1:ntot)+i,:);
        next_x = state(:,i+1);
        if h_2nd >= 0.5
            dt_2nd = 0.1;
            u = repmat(current_u',[ceil(h_2nd/dt_2nd),1]);
            step_record = rk4sim(dynfunc_2nd_burn,current_x',dt_2nd,h_2nd,u);
            ceq(param.nstate*(i-1)+1:param.nstate*i) = next_x - step_record.x(end,:)';
        else
            ceq(param.nstate*(i-1)+1:param.nstate*i) = next_x - singleRK4(dynfunc_2nd_burn,0,h_2nd,current_x,current_u);
        end
    end

    % final condition constraint
    ceq(param.nstate*N+2) = x(N)-param.v_final_ref; % final error in orbital velocity
    ceq(param.nstate*N+3) = x(2*N)-param.fpa_final_ref; % final error in flight path angle
    ceq(param.nstate*N+4) = x(3*N)-param.alt_final_ref; % final error in orbital altitude

    c = nan([N+2,1]);
    % propellant usage constraint
    c(1) = param.m0-state(4,N_1st+1)-param.mp1;
    c(2) = param.m02-state(4,N+1)-param.mp2;
    % above the ground constraint
    c(3:N+2) = -state(3,2:N+1);
 end

%% FUNCTION - Generalized RK4 solver for a dynamic equation
function record = rk4sim(dynFunc,init_cond,dt,tf,u_hist)
    t0 = 0; x = init_cond; forces = [0 0];
    t_hist = t0:dt:tf; x_hist = nan(size(init_cond)); force_hist = zeros([1,2]);
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
% 4 - d, downrange
% 5 - m, total mass

% % u represents control as follows
% 1 - delta TVC angle
% 2 - throttle

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
end
