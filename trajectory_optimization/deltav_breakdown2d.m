clear; clc; close all
addpath("..\parameter_functions")

file_name = "dircol_final1_scaled.mat";
load(file_name,"log_x","log_param");

state_scaling = [1; 1; 1; 1]; time_scaling = 1;
S = log_param.S; maxT_1st = log_param.maxT_1st; maxT_2nd = log_param.maxT_2nd; 
if contains(file_name,"scaled")
    state_scaling = [log_param.scales.speed; 1; log_param.scales.length; log_param.scales.mass];
    time_scaling = log_param.scales.time;
    S = S*log_param.scales.area;
    maxT_1st = maxT_1st*log_param.scales.force;
    maxT_2nd = maxT_2nd*log_param.scales.force;
end

N = log_param.N_1st + log_param.N_2nd;
h_1st = log_x(N*(log_param.nstate+log_param.nctrl)+1)*time_scaling;
h_2nd = log_x(N*(log_param.nstate+log_param.nctrl)+2)*time_scaling;
state = reshape(log_x(1:log_param.nstate*N),[N,log_param.nstate]).*state_scaling';
state_1st = state(1:log_param.N_1st,:);
state_2nd = state(log_param.N_1st+1:N,:);
ctrl = reshape(log_x(log_param.nstate*N+1:N*(log_param.nstate+log_param.nctrl)),[N,log_param.nctrl]);
ctrl_1st = ctrl(1:log_param.N_1st,:);
ctrl_2nd = ctrl(log_param.N_1st+1:N,:);

% coast states - linearly approximated
coast_time = log_x(N*(log_param.nstate+log_param.nctrl)+3);
state_coast = [state_1st(end,:); state_2nd(1,:)];

%% gravity loss
g_1st = g(state_1st(:,3)); g_2nd = g(state_2nd(:,3));
v_gloss_1st = trapz(g_1st.*sin(state_1st(:,2))*h_1st);
v_gloss_2nd = trapz(g_2nd.*sin(state_2nd(:,2))*h_2nd);
v_gloss_coast = trapz(g(state_coast(:,3)).*sin(state_coast(:,2))*coast_time);
v_gloss = v_gloss_1st + v_gloss_coast + v_gloss_2nd;

%% drag loss
rho_1st = rho(state_1st(:,3)); rho_2nd = rho(state_2nd(:,3));
p_1st = P(state_1st(:,3)); p_2nd = P(state_2nd(:,3));
mach_1st = state_1st(:,1)./sqrt(1.4.*p_1st./rho_1st);
for i = 1:length(mach_1st) cd_1st(i) = CD(mach_1st(i),0); end
mach_2nd = state_2nd(:,1)./sqrt(1.4.*p_2nd./rho_2nd);
for i = 1:length(mach_2nd) cd_2nd(i) = CD(mach_2nd(i),0); end
drag_1st = cd_1st'.*0.5.*rho_1st.*state_1st(:,1).^2*S;
drag_2nd = cd_2nd'.*0.5.*rho_2nd.*state_2nd(:,1).^2*S;
v_dloss_1st = trapz(drag_1st./state_1st(:,4)*h_1st);
v_dloss_2nd = trapz(drag_2nd./state_2nd(:,4)*h_2nd);
v_dloss_coast = trapz([drag_1st(end); drag_2nd(1)]/state_2nd(1,4)*h_2nd);
v_dloss = v_dloss_1st + v_dloss_coast + v_dloss_2nd;

%% steering loss
thrust_1st = maxT_1st*ctrl_1st(:,2);
thrust_2nd = maxT_2nd*ctrl_2nd(:,2);
v_sloss_1st = trapz(thrust_1st./state_1st(:,4).*(1-cos(ctrl_1st(:,1)))*h_1st);
v_sloss_2nd = trapz(thrust_2nd./state_2nd(:,4).*(1-cos(ctrl_2nd(:,1)))*h_2nd);
v_sloss = v_sloss_1st + v_sloss_2nd;

%% total delta v
vloss = v_gloss + v_dloss + v_sloss;
deltav_1st = log_param.Isp1*g(0)*log(state_1st(1,4)/state_1st(end,4));
deltav_2nd = log_param.Isp2*g(0)*log(state_2nd(1,4)/state_2nd(end,4));
deltav = deltav_1st + deltav_2nd;

design_dv_1st = log_param.Isp1*g(0)*log(log_param.m0/(log_param.m0-log_param.mp1));
design_dv_2nd = log_param.Isp2*g(0)*log(log_param.m02/(log_param.m02-log_param.mp2));
design_dv = design_dv_1st + design_dv_2nd;
