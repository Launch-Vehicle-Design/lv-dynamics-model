clear; clc; close all

[param, envir] = sysparam();
release_alti = 40000*0.3048;
release_velo = 250.786;

init_cond = [release_velo 0 release_alti 0 param.m0];
% release phase
t_iter = [0,0.1,5]; % 5 second release
[t_rele_record, x_rele_record, x_rele_final, Force_rele_record] = rk4AscentSim(@simplifiedAscent,t_iter,init_cond,1,false,param,envir);
% 1st stage boost phase
t_iter = [t_rele_record(end),0.2,t_rele_record(end)+10]; % arbitrary final time
init_cond_1stb = x_rele_final; init_cond_1stb(2) = 80/180*pi;
[t_1stb_record, x_1stb_record, x_1stb_final, Force_1stb_record] = rk4AscentSim(@simplifiedAscent,t_iter,init_cond_1stb,1,true,param,envir);
% 1st stage coast phase
t_iter = [t_1stb_record(end),0.5,t_1stb_record(end)+350];
init_cond_1stc = x_1stb_final;
[t_1stc_record, x_1stc_record, x_1stc_final, Force_1stc_record] = rk4AscentSim(@simplifiedAscent,t_iter,init_cond_1stc,1,false,param,envir);
% separation
x_2stb_init = x_1stc_final; x_2stb_init(5) = param.m02;
% 2nd stage boost phase
t_iter = [t_1stc_record(end),0.2,t_1stc_record(end)+8]; % arbitrary final time
init_cond_2stb = x_2stb_init; init_cond_2stb(2) = 10/180*pi;
[t_2stb_record, x_2stb_record, x_2stb_final, Force_2stb_record] = rk4AscentSim(@simplifiedAscent,t_iter,init_cond_2stb,2,true,param,envir);
% 2nd stage coast phase
t_iter = [t_2stb_record(end),0.1,t_2stb_record(end)+1];
init_cond_2stc = x_2stb_final;
[t_2stc_record, x_2stc_record, x_2stc_final, Force_2stc_record] = rk4AscentSim(@simplifiedAscent,t_iter,init_cond_2stc,2,false,param,envir);

t_record = [t_rele_record t_1stb_record t_1stc_record t_2stb_record t_2stc_record];
x_record = [x_rele_record; x_1stb_record; x_1stc_record; x_2stb_record; x_2stc_record];
Force_record = [Force_rele_record; Force_1stb_record; Force_1stc_record; Force_2stb_record; Force_2stc_record];

figure; subplot(2,3,1);
plot(t_record, x_record(:,1),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)");
subplot(2,3,2);
plot(t_record, x_record(:,2)*180/pi,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Pitch angle (deg)");
subplot(2,3,3);
plot(t_record, x_record(:,3)/1000,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Altitude (km)");
subplot(2,3,4);
plot(t_record, x_record(:,4)/1000,"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Downrange (km)");
subplot(2,3,5);
plot(t_record, x_record(:,5),"k","LineWidth",1.2); grid on
xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");
subplot(2,3,6); 
plot(x_record(:,4)/1000,x_record(:,3)/1000,"k","LineWidth",1.2); grid on
xlabel("Downrange (km)"); ylabel("Altitude (km)");

figure; plot(t_record,Force_record(:,1),"r","LineWidth",1.5); grid on
hold on; plot(t_record,Force_record(:,2),"k","LineWidth",1.5);
xlabel("Time since Release (s)"); ylabel("Force (N)");
legend("Thrust","Drag")

piece = 50; ealpha = 0:2*pi/piece:2*pi;
ex = envir.Rp*cos(ealpha); ey = envir.Rp*sin(ealpha);
[x,y] = dh2xy(x_record(:,4),x_record(:,3),envir);
figure; plot(x/1000,y/1000,"k","LineWidth",1.5); grid on; hold on
plot(ex/1000,ey/1000,"Color",[0.7 0.7 0.7],"LineWidth",1.2)
xlabel("Earth Centered X (km)"); ylabel("Earth Centered Y (km)");
axis("equal")

sep_time = t_1stc_record(end);
visual2D(t_record,x_record,Force_record,sep_time,param,true,"demo_flight_2.mp4")

function [t,x,x_final,F] = rk4AscentSim(dynFunc,t_iter,init_cond,stage,firing,param,envir)
    t0 = t_iter(1); dt = t_iter(2); tf = t_iter(3);
    if firing
        if stage == 1
            Isp = param.Isp1;
            mp = param.mp1;
            m0 = param.m0;
        elseif stage == 2
            Isp = param.Isp2;
            mp = param.mp2;
            m0 = param.m02;
        end
        tf = t0+Isp;
    end
    tspan = t0:dt:tf; xspan = nan(size(init_cond)); Fspan = nan([1,2]);
    x = init_cond;
    for i = 1:size(tspan,2)
        t = tspan(i);
        xspan(i,:) = x;
        [k1,F1] = dynFunc(t,x,stage,firing,param,envir);
        [k2,F2] = dynFunc(t+dt/2,x+k1*dt/2,stage,firing,param,envir);
        [k3,F3] = dynFunc(t+dt/2,x+k2*dt/2,stage,firing,param,envir);
        [k4,F4] = dynFunc(t+dt,x+k3*dt,stage,firing,param,envir);
        x = x + 1/6*dt*(k1+2*k2+2*k3+k4);
        Fspan(i,:) = 1/6*(F1+2*F2+2*F3+F4);
        if firing && x(5) < m0-0.99*mp
            break
        end
    end

    t = tspan(1:i); x_final = x; x = xspan; F = Fspan;
end

function [x,y] = dh2xy(downrange, altitude, envir)
    turn_angle = downrange./envir.Rp;
    x = (envir.Rp+altitude).*sin(turn_angle);
    y = (envir.Rp+altitude).*cos(turn_angle);
end