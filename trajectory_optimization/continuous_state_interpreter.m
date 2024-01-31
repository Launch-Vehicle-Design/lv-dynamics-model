function state = continuous_state_interpreter(t,log_x,log_param)
% t, post drop phase time, zero at the end of the drop
% file_name (optional), solution mat file name
% author: Chenming Fan

x0 = [0.032136009063672; 1.396263401595464; 0.001898721754591; 1];

state = nan([4,1]);

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
end