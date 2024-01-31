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

% environmental parameter
param.mu = 3.986004418e14;
param.earthR = 6357e3;
param.omega_len = 2*pi/(24*3600);
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

param.TtoW_1st = 1.6; param.TtoW_2nd = 0.7;
param.maxT_1st = param.TtoW_1st*param.m0*9.80665;
param.maxT_2nd = param.TtoW_2nd*param.m02*9.80665;

init_cond = [release_velo 0 release_alti 0 param.m0];
% PHASE 1 - uncontrolled drop phase
drop_dt = 0.01; drop_time = 5;
drop_dynfun = @(t,x,u) rocket(t,x,u,param,0);
drop_record = rk4sim(drop_dynfun,init_cond,drop_dt,drop_time,[0 0]);

% % PHASE 1 - plotter
% figure; subplot(2,3,1); plot(drop_record.t, drop_record.x(:,1),"k","LineWidth",1.2); grid on
% xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
% subplot(2,3,2); plot(drop_record.t, drop_record.x(:,2)*180/pi,"k","LineWidth",1.2); grid on
% xlabel("Time since Release (s)"); ylabel("Pitch angle (deg)");
% subplot(2,3,3); plot(drop_record.t, drop_record.x(:,3)/1000,"k","LineWidth",1.2); grid on
% xlabel("Time since Release (s)"); ylabel("Altitude (km)");
% subplot(2,3,4); plot(drop_record.t, drop_record.x(:,4)/1000,"k","LineWidth",1.2); grid on
% xlabel("Time since Release (s)"); ylabel("Downrange (km)");
% subplot(2,3,5); plot(drop_record.t, drop_record.x(:,5),"k","LineWidth",1.2); grid on
% xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");

%% Trajectory Optimization - direct trajectory optimization
v_final_ref = 7754.1; alt_final_ref = 251460; fpa_final_ref = 0;
% xopt represents state and control trajectory of the system
% with total length of 7*N+1
% 1:5N - state from x11...xN1:x51...x5N
% 5N+1:7N - control from u10...u1N-1,u20...u2N-1
% 7N+1:7N+2 - the step length for both 1st and 2nd stage
% 7N+3 - the coast time at stage separation

% cost function weighting
mp_weight = 1e-5; fpaf_weight = 1000; vf_weight = 1e-6; altf_weight = 1e-8;
dtvc_weight_1st = 100; dthrotl_weight_1st = 1;
dtvc_weight_2nd = 1000; dthrotl_weight_2nd = 1;

% design control bounds
tvc_bounds = [-15 15]/180*pi;
throtl_bounds = [0.4 1.01];

N_1st = 200; N_2nd = 250;
N = N_1st+N_2nd; x0 = drop_record.x(end,:)';
x0(2) = 80*pi/180;

% initial guess propagation
init_x = zeros([7*N+3,1]); 
init_x(7*N+1) = 0.85; 
init_x(7*N+2) = 1.2; 
init_x(7*N+3) = 5;
param.init_ind_ptr = 0:N:7*N;
current_x = x0; init_u = [0 0.8];
% init_x(5*N+1:7*N) = reshape(repmat(init_u,[N,1]),[2*N,1]);
for i = 1:N_1st
    u = init_u;
    u = [max(min(-1*(current_x(2)-50/180*pi),tvc_bounds(2)*0.99),tvc_bounds(1)*0.99) init_u(2)];
    if param.m0 - current_x(5) >= param.mp1 u(2) = 0; end
    dynfunc_1st = @(t,x,u) rocket(t,x,u,param,1);
    current_x = singleRK4(dynfunc_1st,0,init_x(7*N+1),current_x,u);
    init_x(param.init_ind_ptr(1:5)+i) = current_x;
    init_x(param.init_ind_ptr(6:7)+i) = u;
end
current_x(5) = param.m02;
dynfunc_2nd = @(t,x,u) rocket(t,x,u,param,2);
stg_sep_record = rk4sim(dynfunc_2nd,current_x',0.1,init_x(7*N+3),[0 0]);
current_x = stg_sep_record.x(end,:)';
init_u = [0 1];
init_x(param.init_ind_ptr(1:5)+N_1st+1) = current_x;
init_x(param.init_ind_ptr(6:7)+N_1st+1) = init_u;
for i = N_1st+2:N_1st+N_2nd
    u = init_u;
    u = [max(min(-1*(current_x(2)-20/180*pi),tvc_bounds(2)*0.99),tvc_bounds(1)*0.99) init_u(2)];
    if param.m02 - current_x(5) >= param.mp2 u(2) = 0; end
    current_x = singleRK4(dynfunc_2nd,0,init_x(7*N+2),current_x,u);
    init_x(param.init_ind_ptr(1:5)+i) = current_x;
    init_x(param.init_ind_ptr(6:7)+i) = u;
end
plotter(init_x,drop_time,N_1st,N_2nd,param);

% cost function
cost_func = @(xopt) -mp_weight*xopt(5*N)^2 + ... % minimize fuel used
    vf_weight*(xopt(N)-v_final_ref)^2 + ... % final error in orbital velocity
    fpaf_weight*(xopt(2*N)-fpa_final_ref)^2 + ... % final error in flight path angle
    altf_weight*(xopt(3*N)-alt_final_ref)^2 + ... % final error in orbital altitude
    xopt(7*N+1)*(dtvc_weight_1st*sum(diff(xopt(5*N+1:5*N+1+N_1st)).^2) + dthrotl_weight_1st*sum(diff(xopt(6*N+1:6*N+N_1st)).^2)) + ...
    xopt(7*N+2)*(dtvc_weight_2nd*sum(diff(xopt(5*N+1+N_1st:6*N)).^2) + dthrotl_weight_2nd*sum(diff(xopt(6*N+1+N_1st:7*N)).^2)); % minimize control frequency

cost_func(init_x)

% state and control xopt bounds
low_bound = [-inf*ones([5*N,1]);
    tvc_bounds(1)*ones([N,1]); throtl_bounds(1)*ones([N,1]);
    0; 0; 0];
upp_bound = [inf*ones([5*N,1]);
    tvc_bounds(2)*ones([N,1]); throtl_bounds(2)*ones([N,1]);
    5; 5; 100];

% constraint function
nonl_con_func = @(x) nonl_con_sep(x,N_1st,N_2nd,x0,param);

% optimization solver setup
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
    "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',1e8,'MaxIterations',1e7,...
    EnableFeasibilityMode=true);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
%     "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',1e8,'MaxIterations',1e7);
x = fmincon(cost_func,init_x,[],[],[],[],low_bound,upp_bound,nonl_con_func,options);
plotter(x,drop_time,N_1st,N_2nd,param);

%% FUNCTION - State and Control Plotter
function plotter(x,drop_time,N_1st,N_2nd,param)
N = N_1st + N_2nd;
t_1st = drop_time:x(7*N+1):drop_time+x(7*N+1)*(N_1st-1);
t_2nd = drop_time+x(7*N+1)*N_1st:x(7*N+2):drop_time+x(7*N+1)*N_1st+x(7*N+2)*(N_2nd-1);
figure; subplot(2,3,1); plot([t_1st t_2nd], x(param.init_ind_ptr(1)+(1:N)),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
subplot(2,3,2); plot([t_1st t_2nd], x(param.init_ind_ptr(2)+(1:N))*180/pi,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Pitch angle (deg)");
subplot(2,3,3); plot([t_1st t_2nd], x(param.init_ind_ptr(3)+(1:N))/1000,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Altitude (km)");
subplot(2,3,4); plot([t_1st t_2nd], x(param.init_ind_ptr(4)+(1:N))/1000,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Downrange (km)");
subplot(2,3,5); plot([t_1st t_2nd], x(param.init_ind_ptr(5)+(1:N)),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");
subplot(2,3,6); plot([t_1st t_2nd], x(param.init_ind_ptr(6)+(1:N))*180/pi,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Gimbal Angle (deg)");
hold on; plot([t_1st t_2nd], x(param.init_ind_ptr(7)+(1:N)),"k","LineWidth",1.2);
end

%% FUNCTION - Constraint function
% function [c,ceq] = nonl_con(x,N,x0,param)
%     h = x(end); state = [x0 reshape(x(1:5*N),[5,N])];
%     ceq = nan([5*N,1]);
%     init_ind_ptr = 1:N:7*N+1;
%     stage = 1;
%     for i = 1:N
%         current_x = state(:,i);
%         current_u = x(init_ind_ptr(6:7)+i,:);
%         next_x = state(:,i+1);
%         % decide the stage seperation
%         if param.m0 - (current_x(5) + next_x(5))/2 <= param.mp1
%             stage = 2;
%         end
%         dynfunc = @(t,x,u) rocket(t,x,u,param,stage);
%         ceq(5*(i-1)+1:5*i) = next_x - singleRK4(dynfunc,0,h,current_x,current_u);
%     end
%     c = [];
% end

function [c,ceq] = nonl_con_sep(x,N_1st,N_2nd,x0,param)
    N = N_1st + N_2nd;
    h_1st = x(7*N+1); h_2nd = x(7*N+2);
    coast_t = x(7*N+3);

    state = [x0 reshape(x(1:5*N),[N,5])'];
    ceq = nan([5*N+2,1]);
    init_ind_ptr = param.init_ind_ptr;

    % equality constraint from dynamics of the 1st stage burn
    % from x0 to x_N1st and u0 to u_N1st-1
    % total of N1st+1 states and N1st controls
    for i = 1:N_1st
        current_x = state(:,i);
        current_u = x(init_ind_ptr(6:7)+i,:);
        next_x = state(:,i+1);
        dynfunc_1st_burn = @(t,x,u) rocket(t,x,u,param,1);
        ceq(5*(i-1)+1:5*i) = next_x - singleRK4(dynfunc_1st_burn,0,h_1st,current_x,current_u);
    end

    % stage separation constraints from x_N1st to x_N1st+1 reset vehicle mass
    coast_dt = 0.1;
    stage_sep_x = [next_x(1:4); param.m02];
    dynfunc_stg_sep = @(t,x,u) rocket(t,x,u,param,2);
    stg_sep_record = rk4sim(dynfunc_stg_sep,stage_sep_x',coast_dt,coast_t,[0 0]);
    ceq(5*N_1st+1:5*(N_1st+1)) = state(:,N_1st+2) - stg_sep_record.x(end,:)';
    % ensure initial zero TVC actuation for 2nd stage
    ceq(5*N+1:5*N+2) = x(init_ind_ptr(6:7)+N_1st+1,:);
    
    % equality constraint from dynamics of the 2nd stage burn
    % from x_N1st+1 to x_N and uN1st+1 to u_N-1
    % total of N2nd states and N2nd-1 controls
    for i = N_1st+2:N_1st+N_2nd
        current_x = state(:,i);
        current_u = x(init_ind_ptr(6:7)+i,:);
        next_x = state(:,i+1);
        dynfunc_2nd_burn = @(t,x,u) rocket(t,x,u,param,2);
        ceq(5*(i-1)+1:5*i) = next_x - singleRK4(dynfunc_2nd_burn,0,h_2nd,current_x,current_u);
        % ceq(5*(i-1)+1:5*i) = next_x - h_2nd*dynfunc_2nd_burn(0,current_x,current_u);
    end

    % propellant usage constraint
    c(1) = param.m0-state(5,N_1st+1)-param.mp1;
    c(2) = param.m02-state(5,N+1)-param.mp2;
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

v = x(1); fpa = x(2); h = x(3); m = x(5);
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
dxdt(4) = v*cos(fpa)/(1+h/Rp);
dxdt(5) = -mdot;
forces = [T D];
end
