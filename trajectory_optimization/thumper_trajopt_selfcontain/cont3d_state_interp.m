function state = cont3d_state_interp(t,log_x,log_param)
% t, post drop phase time, zero at the end of the drop
% file_name (optional), solution mat file name
% author: Chenming Fan

% x0 = [0.032136009063672; 1.396263401595464; 0.001898721754591; 1];
x0 = [0.173977866930877; -0.000196452384873043; 0.986677513325079;
    0.00555850916780049; -0.000636441747248759; 0.0315238679430946; 1];

state = nan([log_param.nstate,1]);

dynfunc_1st_burn = @(t,x,u) rocket(t,x,u,log_param,1);
dynfunc_2nd_burn = @(t,x,u) rocket(t,x,u,log_param,2);
state_scaling = [1; 1; 1; 1]; time_scaling = 1;
if log_param.scales.length ~= 1
    state_scaling = [log_param.scales.speed; 1; log_param.scales.length; log_param.scales.mass];
    time_scaling = log_param.scales.time;
end
t = t/time_scaling;
log_N1st = log_param.N_1st; log_N2nd = log_param.N_2nd;
log_N = log_N1st + log_N2nd; log_ntot = log_param.nstate + log_param.nctrl;
log_h1st = log_x(log_ntot*log_N+1);
log_h2nd = log_x(log_ntot*log_N+2);
log_tcoast = log_x(log_ntot*log_N+3);
log_state = [x0 reshape(log_x(1:log_param.nstate*log_N),[log_N,log_param.nstate])'];
log_control = [[0;1] reshape(log_x(log_param.nstate*log_N+1:log_ntot*log_N),[log_N,log_param.nctrl])'];

% Hermite spline fit
C = @(h,x) [1 0 0 0; 0 1 0 0; -3/h^2 -2/h 3/h^2 -1/h; 2/h^3 1/h^2 -2/h^3 1/h^2]*x;
xspline = @(t,C,h) C(1,floor(t/h)+1)+C(2,floor(t/h)+1).*mod(t,h)+...
    C(3,floor(t/h)+1).*mod(t,h).^2+C(4,floor(t/h)+1).*mod(t,h).^3;
if t < log_h1st*log_N1st
    % 1st stage phase
    ind_t = ceil(t/log_h1st);
    residue_t = mod(t,log_h1st);
    % derivatives at end states
    log_xdot_1st = dynfunc_1st_burn(0,log_state(:,ind_t:ind_t+1),log_control(:,ind_t:ind_t+1));
    for state_ind = 1:log_param.nstate
        xk = log_state(state_ind,ind_t); xkp1 = log_state(state_ind,ind_t+1);
        xdotk = log_xdot_1st(state_ind,1); xdotkp1 = log_xdot_1st(state_ind,2);
        C_1st = C(log_h1st,[xk; xdotk; xkp1; xdotkp1]);
        state(state_ind) = xspline(residue_t,C_1st,log_h1st)';
    end

elseif t <= log_h1st*log_N1st + log_tcoast
    % coast phase
    residue_t = t-log_h1st*log_N1st; dt = 0.01;
    init_cond = [log_state(1:log_param.nstate-1,1+log_N1st); log_state(log_param.nstate,log_N1st+2)];
    record = rk4sim(dynfunc_2nd_burn,init_cond,dt/log_param.scales.time,residue_t,[0 0]);
    state = record.x(:,end);

elseif t < log_h1st*log_N1st + log_tcoast + log_h2nd*(log_N2nd-1)
    % 2nd stage phase
    ind_t = ceil((t-log_h1st*log_N1st-log_tcoast)/log_h2nd);
    residue_t = mod(t-log_h1st*log_N1st-log_tcoast,log_h2nd);
    % derivatives at end states
    log_xdot_2nd = dynfunc_2nd_burn(0,log_state(:,log_N1st+ind_t+1:log_N1st+ind_t+2),log_control(:,log_N1st+ind_t+1:log_N1st+ind_t+2));
    for state_ind = 1:log_param.nstate
        xk = log_state(state_ind,log_N1st+ind_t+1); xkp1 = log_state(state_ind,log_N1st+2+ind_t);
        xdotk = log_xdot_2nd(state_ind,1); xdotkp1 = log_xdot_2nd(state_ind,2);
        C_2nd = C(log_h2nd,[xk; xdotk; xkp1; xdotkp1]);
        state(state_ind) = xspline(residue_t,C_2nd,log_h2nd)';
    end
else
    state = log_state(:,end);
end
state = state.*state_scaling;

%% FUNCTION - Generalized RK4 solver for a dynamic equation
    function record = rk4sim(dynFunc,init_cond,dt,tf,u_hist)
        t0 = 0; x = init_cond; forces = [0 0];
        t_hist = t0:dt:tf; x_hist = init_cond; force_hist = zeros([1,2]);
        for i = 1:size(t_hist,2)
            tnow = t_hist(i); x_hist(:,i) = x; force_hist(i,:) = forces;
            if i <= size(u_hist,1) u = u_hist(i,:); else u = [0 0]; end
            [x,forces] = singleRK4(dynFunc,tnow,dt,x,u);
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
end