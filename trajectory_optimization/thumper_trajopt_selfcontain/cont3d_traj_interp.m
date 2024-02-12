function cont_traj = cont3d_traj_interp(t,log_x,log_param)
% t, post drop phase time, zero at the end of the drop
% file_name (optional), solution mat file name
% author: Chenming Fan

% x0 = [0.173977866930877; -0.000196452384873043; 0.986677513325079;
%     0.00555850916780049; -0.000636441747248759; 0.0315238679430946; 1];
x0 = log_param.x0;

state = nan([log_param.nstate,1]);
ctrl = nan([log_param.nctrl,1]);
zero_ctrl = zeros([log_param.nctrl,1]);

dynfunc_1st_burn = @(t,x,u) rocket(t,x,u,log_param,1);
dynfunc_2nd_burn = @(t,x,u) rocket(t,x,u,log_param,2);

state_scaling = ones([log_param.nstate,1]); time_scaling = 1;
if log_param.earthR == 1
    state_scaling = [log_param.scales.length*ones([3 1]);
        log_param.scales.speed*ones([3 1]);
        log_param.scales.mass];
    time_scaling = log_param.scales.time;
end
t = t/time_scaling;

log_t_drop = log_param.drop_time;
log_N1st = log_param.N_1st; log_N2nd = log_param.N_2nd;
log_N = log_N1st + log_N2nd; log_ntot = log_param.nstate + log_param.nctrl;
log_h1st = log_x(log_ntot*log_N+1);
log_h2nd = log_x(log_ntot*log_N+2);
log_tcoast = log_x(log_ntot*log_N+3);
log_state = [x0 reshape(log_x(1:log_param.nstate*log_N),[log_N,log_param.nstate])'];
log_control = [[x0(4:6)/norm(4:6);1] reshape(log_x(log_param.nstate*log_N+1:log_ntot*log_N),[log_N,log_param.nctrl])'];

% Hermite spline fit
C = @(h,x) [1 0 0 0; 0 1 0 0; -3/h^2 -2/h 3/h^2 -1/h; 2/h^3 1/h^2 -2/h^3 1/h^2]*x;
xspline = @(t,C,h) C(1,floor(t/h)+1)+C(2,floor(t/h)+1).*mod(t,h)+...
    C(3,floor(t/h)+1).*mod(t,h).^2+C(4,floor(t/h)+1).*mod(t,h).^3;
ulinear = @(t,uk,ukp1,h) mod(t,h).*(ukp1(floor(t/h)+1)-uk(floor(t/h)+1))./h + uk(floor(t/h)+1);

% mission time line
mission_tline = [log_t_drop;
    log_t_drop + log_h1st*log_N1st;
    log_t_drop + log_h1st*log_N1st + log_tcoast;
    log_t_drop + log_h1st*log_N1st + log_tcoast + log_h2nd*(log_N2nd-1)];

% drop phase
if t <= mission_tline(1)
    residue_t = t; dt = 0.01;
    init_cond = log_param.init_cond;
    record = rk4sim(dynfunc_1st_burn,init_cond,dt/time_scaling,residue_t,zero_ctrl);
    state = record.x(:,end);
    ctrl = [state(4:6)/norm(state(4:6)); 0];

% 1st stage phase
elseif t < mission_tline(2)
    ind_t = ceil((t-mission_tline(1))/log_h1st);
    residue_t = mod(t-mission_tline(1),log_h1st);
    % derivatives at end states
    log_xdot_1st = dynfunc_1st_burn(0,log_state(:,ind_t:ind_t+1),log_control(:,ind_t:ind_t+1));
    for state_ind = 1:log_param.nstate
        xk = log_state(state_ind,ind_t); xkp1 = log_state(state_ind,ind_t+1);
        xdotk = log_xdot_1st(state_ind,1); xdotkp1 = log_xdot_1st(state_ind,2);
        C_1st = C(log_h1st,[xk; xdotk; xkp1; xdotkp1]);
        state(state_ind) = xspline(residue_t,C_1st,log_h1st)';
    end
    for ctrl_ind = 1:log_param.nctrl
        uk = log_control(ctrl_ind,ind_t); ukp1 = log_control(ctrl_ind,ind_t+1);
        ctrl(ctrl_ind) = ulinear(residue_t,uk,ukp1,log_h1st)';
    end

% coast phase
elseif t <= mission_tline(3) 
    residue_t = t-mission_tline(2); dt = 0.01;
    init_cond = [log_state(1:log_param.nstate-1,1+log_N1st); log_state(log_param.nstate,log_N1st+2)];
    record = rk4sim(dynfunc_2nd_burn,init_cond,dt/log_param.scales.time,residue_t,zero_ctrl);
    state = record.x(:,end);
    ctrl = [state(4:6)/norm(state(4:6)); 0];

% 2nd stage phase
elseif t < mission_tline(4)
    ind_t = ceil((t-mission_tline(3))/log_h2nd);
    residue_t = mod(t-mission_tline(3),log_h2nd);
    % derivatives at end states
    log_xdot_2nd = dynfunc_2nd_burn(0,log_state(:,log_N1st+ind_t+1:log_N1st+ind_t+2),log_control(:,log_N1st+ind_t+1:log_N1st+ind_t+2));
    for state_ind = 1:log_param.nstate
        xk = log_state(state_ind,log_N1st+ind_t+1); xkp1 = log_state(state_ind,log_N1st+2+ind_t);
        xdotk = log_xdot_2nd(state_ind,1); xdotkp1 = log_xdot_2nd(state_ind,2);
        C_2nd = C(log_h2nd,[xk; xdotk; xkp1; xdotkp1]);
        state(state_ind) = xspline(residue_t,C_2nd,log_h2nd)';
    end
    for ctrl_ind = 1:log_param.nctrl
        uk = log_control(ctrl_ind,log_N1st+ind_t); ukp1 = log_control(ctrl_ind,log_N1st+2+ind_t);
        ctrl(ctrl_ind) = ulinear(residue_t,uk,ukp1,log_h2nd)';
    end
else
    state = log_state(:,end);
end
state = state.*state_scaling;

cont_traj.state = state;
cont_traj.ctrl = ctrl;
end