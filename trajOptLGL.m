clear; clc; close all
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

%% Range Band Plot
% initState = [[5e5; 1e5; 5e4]*0.3048;... reference
%             [-3000; 0; 0]*0.3048];

insert_downRange = -7.8e5:1e5:-7.7e5;
crossRange = 0;
release_altitude = 12.192e3;
altitude = 251.46e3;
targ_velocity = [0; 7754.1; 0];
init_velocity = [250; 0; 0];

% constant setup
% evaluating condition
cond.mu = 3.986004418e14;  
cond.R = 6378e3;
cond.dt = 0.1;
cond.targetLocation = [cond.R+altitude; 0; 0];
cond.tolerance = 1;
cond.Gc = -cond.mu*cond.targetLocation/norm(cond.targetLocation)^3;
cond.Gi = -cond.mu/norm(cond.targetLocation)^3*...
    (eye(3)-3/norm(cond.targetLocation)^2*cond.targetLocation.^2.*eye(3));
cond.tgoGuess = 10;
cond.gamma = 1500;
cond.options = optimset('TolX',1e-4);
cond.v0 = targ_velocity;
cond.vf = [7754.1; 0; 0];
cond.tout = 5000;

min_dv = 1e9;
for i = 1:size(insert_downRange,2)
    initState = [-insert_downRange(i); crossRange; release_altitude; init_velocity];
    dsouzaInitState = [-insert_downRange(i); crossRange; release_altitude-altitude; init_velocity-cond.vf];
    [history_dsouza, perIndex_dsouza] = dsouzaGuidLaw(dsouzaInitState,cond);
    % [history_dsouza, perIndex_dsouza] = dsouzaOrbitInsertGuidLaw(dsouzaInitState,cond);
    cond.tgoGuess = max(history_dsouza.time);
    % [history_linear, perIndex_linear] = linearOrbitInsertGuidLaw(initState,cond);

    sanitized_a = history_dsouza.a;
    sanitized_a(isnan(sanitized_a)) = 0;
    isp = 300; g = 9.80665;
    m0 = 1; mf = 1;
    for ind = 1:length(history_dsouza.time)
        mdot = norm(sanitized_a(ind,:))*mf/(isp*g);
        mf = mf - mdot*cond.dt;
        m(ind) = mf;
    end
    dv = isp*g*log(1./m)
    [val,sep_ind] = min(abs(dv-4666));
    max(history_dsouza.time)

    if dv < min_dv
        min_dv = dv;
    end

    plotHistory(history_dsouza);
    % plotHistory(history_linear);

    figure;
    grd_x = @(x) x+7754.1*history_dsouza.time+insert_downRange(i);
    ground_x = grd_x(history_dsouza.x(:,1));
    plot(ground_x/0.3048,(history_dsouza.x(:,3)+altitude)/0.3048,"LineWidth",1.2); hold on;
    yline(100e3/0.3048,"--");
    scatter(ground_x(sep_ind)/0.3048,(history_dsouza.x(sep_ind,3)+altitude)/0.3048,"LineWidth",1.2)
    text(100e3,825000,0,"Desired Orbit Altitude at 825,000 ft")
    xlabel("downrange (ft)"); ylabel("altitude (ft)"); axis("equal"); grid on;
    legend("Trajectory", "Karman Line", "1st Step Separation","interpreter","latex");
end

function plotHistory(history)
    figure;
    % acceleration plot location 1,4,7
    subplot(3,3,1); plot(history.time, history.a(:,1)); grid on
    title("ax"); xlabel("time (s)"); ylabel("x acce ($m/s^2$)");
    subplot(3,3,4); plot(history.time, history.a(:,2)); grid on
    title("ay"); xlabel("time (s)"); ylabel("y acce ($m/s^2$)");
    subplot(3,3,7); plot(history.time, history.a(:,3)); grid on
    title("az"); xlabel("time (s)"); ylabel("z acce ($m/s^2$)");
    
    % velocity plot location 2,5,8
    subplot(3,3,2); plot(history.time, history.v(:,1)); grid on
    title("vx"); xlabel("time (s)"); ylabel("x velo (m/s)");
    subplot(3,3,5); plot(history.time, history.v(:,2)); grid on
    title("vy"); xlabel("time (s)"); ylabel("y velo (m/s)");
    subplot(3,3,8); plot(history.time, history.v(:,3)); grid on
    title("vz"); xlabel("time (s)"); ylabel("z velo (m/s)");
    
    % position plot location 3,6,9
    subplot(3,3,3); plot(history.time, history.x(:,1)); grid on
    title("x"); xlabel("time (s)"); ylabel("x (m)");
    subplot(3,3,6); plot(history.time, history.x(:,2)); grid on
    title("y"); xlabel("time (s)"); ylabel("y (m)");
    subplot(3,3,9); plot(history.time, history.x(:,3)); grid on
    title("z"); xlabel("time (s)"); ylabel("z (m)");
end
