file_name = "thumper_straj_cg3dof.mat";
load(file_name,"log_x","log_param");

state_scaling = ones([7 1]); time_scaling = 1;
S = log_param.S; maxT_1st = log_param.maxT_1st; maxT_2nd = log_param.maxT_2nd;
Re = log_param.earthR;
if contains(file_name,"scaled") || contains(file_name,"straj")
    state_scaling = [
        log_param.scales.length*ones([3 1]);
        log_param.scales.speed*ones([3 1]);
        log_param.scales.mass
        ];
    time_scaling = log_param.scales.time;
    S = S*log_param.scales.area;
    maxT_1st = maxT_1st*log_param.scales.force;
    maxT_2nd = maxT_2nd*log_param.scales.force;
    Re = Re*log_param.scales.length;
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
e_r1st = state_1st(:,1:3)'./vecnorm(state_1st(:,1:3)');
e_r2nd = state_2nd(:,1:3)'./vecnorm(state_2nd(:,1:3)');
e_v1st = state_1st(:,4:6)'./vecnorm(state_1st(:,4:6)');
e_v2nd = state_2nd(:,4:6)'./vecnorm(state_2nd(:,4:6)');
e_u1st = ctrl_1st(:,1:3)'./vecnorm(ctrl_1st(:,1:3)');
e_u2nd = ctrl_2nd(:,1:3)'./vecnorm(ctrl_2nd(:,1:3)');
e_rcst = state_coast(:,1:3)'./vecnorm(state_coast(:,1:3)');
e_vcst = state_coast(:,4:6)'./vecnorm(state_coast(:,4:6)');
alt_1st = vecnorm(state_1st(:,1:3)')-Re; alt_2nd = vecnorm(state_2nd(:,1:3)')-Re;

%% gravity loss
mu = log_param.mu*log_param.scales.gravparam;
g_1st = mu./vecnorm(state_1st(:,1:3)').^2;
g_2nd = mu./vecnorm(state_2nd(:,1:3)').^2;
v_gloss_1st = trapz(g_1st.*dot(e_r1st,e_v1st)*h_1st);
v_gloss_2nd = trapz(g_2nd.*dot(e_r2nd,e_v2nd)*h_2nd);
v_gloss_coast = trapz(g(vecnorm(state_coast(:,1:3)')-Re).*dot(e_rcst,e_vcst)*coast_time);
v_gloss = v_gloss_1st + v_gloss_coast + v_gloss_2nd;

%% drag loss
atmo_profile_1st = atmo(alt_1st); atmo_profile_2nd = atmo(alt_2nd);
rho_1st = atmo_profile_1st.rho; rho_2nd = atmo_profile_2nd.rho;
p_1st = atmo_profile_1st.P; p_2nd = atmo_profile_2nd.P;
% mach_1st = vecnorm(state_1st(:,4:6)')./sqrt(1.4.*p_1st./rho_1st);
% for i = 1:length(mach_1st) cd_1st(i) = CD(mach_1st(i),0); end
% mach_2nd = vecnorm(state_2nd(:,4:6)')./sqrt(1.4.*p_2nd./rho_2nd);
% for i = 1:length(mach_2nd) cd_2nd(i) = CD(mach_2nd(i),0); end
cd_1st = 0.5; cd_2nd = 0.5;
drag_1st = cd_1st.*0.5.*rho_1st.*vecnorm(state_1st(:,4:6)').^2*S;
drag_2nd = cd_2nd.*0.5.*rho_2nd.*vecnorm(state_2nd(:,4:6)').^2*S;
v_dloss_1st = trapz(drag_1st./state_1st(:,7)'*h_1st);
v_dloss_2nd = trapz(drag_2nd./state_2nd(:,7)'*h_2nd);
v_dloss_coast = trapz([drag_1st(end); drag_2nd(1)]/state_2nd(1,7)*h_2nd);
v_dloss = v_dloss_1st + v_dloss_coast + v_dloss_2nd;

%% steering loss
thrust_1st = maxT_1st*ctrl_1st(:,4)';
thrust_2nd = maxT_2nd*ctrl_2nd(:,4)';
v_sloss_1st = trapz(thrust_1st./state_1st(:,7)'.*(1-dot(e_u1st,e_v1st))*h_1st);
v_sloss_2nd = trapz(thrust_2nd./state_2nd(:,7)'.*(1-dot(e_u2nd,e_v2nd))*h_2nd);
v_sloss = v_sloss_1st + v_sloss_2nd;

%% total delta v
vloss = v_gloss + v_dloss + v_sloss;
deltav_1st = log_param.Isp1*g(0)*log(state_1st(1,4)/state_1st(end,4));
deltav_2nd = log_param.Isp2*g(0)*log(state_2nd(1,4)/state_2nd(end,4));
deltav = deltav_1st + deltav_2nd;

design_dv_1st = log_param.Isp1*g(0)*log(log_param.m0/(log_param.m0-log_param.mp1));
design_dv_2nd = log_param.Isp2*g(0)*log(log_param.m02/(log_param.m02-log_param.mp2));
design_dv = design_dv_1st + design_dv_2nd;
