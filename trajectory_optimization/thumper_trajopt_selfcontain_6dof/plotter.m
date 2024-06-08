%% FUNCTION - State and Control Plotter
function axes = plotter(x,param)
    N_1st = param.N_1st; N_2nd = param.N_2nd;
    N = N_1st + N_2nd;
    ntot = param.nstate + param.nctrl;
    t_1st = 0:x(ntot*N+1):x(ntot*N+1)*(N_1st-1);
    t_2nd = x(ntot*N+1)*N_1st:x(ntot*N+2):x(ntot*N+1)*N_1st+x(ntot*N+2)*(N_2nd-1);
    unsc_t = [t_1st t_2nd]*param.scales.time;

    states = reshape(x(1:param.nstate*N),[N param.nstate])'.*param.scales.state;
    ctrls = reshape(x(param.nstate*N+1:N*ntot),[N param.nctrl])';
    % quatornion separation
    q0 = states(7,:); q1 = states(8,:); q2 = states(9,:); q3 = states(10,:);
    rn = vecnorm(states(1:3,:)); vn = vecnorm(states(4:6,:));
    fpa = pi/2-acos(dot(states(1:3,:),states(4:6,:))./rn./vn);
    % convert quaternion to rotation matrix
    % C(1,1) = q0*q0+q1*q1-q2*q2-q3*q3;
    % C(1,2) = 2*(q1*q2+q0*q3);
    % C(1,3) = 2*(q1*q3-q0*q2);
    % C(2,1) = 2*(q1*q2-q0*q3);
    % C(2,2) = q0*q0-q1*q1+q2*q2-q3*q3;
    % C(2,3) = 2*(q2*q3+q0*q1);
    % C(3,1) = 2*(q1*q3 + q0*q2);
    % C(3,2) = 2*(q2*q3-q0*q1);
    % C(3,3) = q0*q0-q1*q1-q2*q2+q3*q3;
    xiner = [q0.^2 + q1.^2 - q2.^2 - q3.^2;
        2*(q1.*q2 - q0.*q3); 2*(q1.*q3 + q0.*q2)];
    aoa = acos(dot(xiner,states(4:6,:))./vn);

    show_num = min(30,N); show_ind = 1:floor(N/show_num):N;
    show_t = unsc_t(show_ind);
    
    subplot_size = [4 4]; subplot_ind = 1;
    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind);
    plot3(states(1,show_ind),states(2,show_ind),states(3,show_ind),"k","LineWidth",2); hold on;
    quiver3(states(1,show_ind),states(2,show_ind),states(3,show_ind),xiner(1,show_ind),xiner(2,show_ind),xiner(3,show_ind));
    xl = xlim; yl = ylim; zl = zlim; [X,Y,Z] = sphere(100);
    surf(X*param.scales.length,Y*param.scales.length,Z*param.scales.length,"EdgeColor","none","FaceColor",[0.7 0.7 0.7]);
    axis equal; xlim(xl); ylim(yl); zlim(zl);
    xlabel("x"); ylabel("y"); zlabel("z"); hold off
    subplot_ind = subplot_ind + 1;

    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind);
    plot(show_t,rn(show_ind)-param.earthR*param.scales.length,"k","LineWidth",1.2); grid on;
    xlabel("Time since Release (s)"); ylabel("Altitude (m)"); xlim([0,600]); hold off;
    subplot_ind = subplot_ind + 1;

    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind);
    plot(show_t,vn(show_ind),"k","LineWidth",1.2); grid on;
    xlabel("Time since Release (s)"); ylabel("Speed (m/s)"); xlim([0,600]); hold off;
    subplot_ind = subplot_ind + 1;

    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind);
    plot(show_t, states(14,show_ind),"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("Mass (kg)");
    xlim([0,600]); ylim([0,1850]);
    subplot_ind = subplot_ind + 1;

    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind);
    plot(show_t, states(15,show_ind),"k","LineWidth",1.2); hold on; grid on
    plot(show_t, states(16,show_ind),"k","LineWidth",1.2);
    xlabel("Time since Release (s)"); ylabel("LiqProp Mass (kg)"); 
    xlim([0,600]); hold off;
    subplot_ind = subplot_ind + 1;

    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind);
    plot(show_t,fpa(show_ind)*180/pi,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("FPA (deg)");
    xlim([0,600]);
    subplot_ind = subplot_ind + 1;

    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind);
    plot(show_t,aoa(show_ind)*180/pi,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("AoA (deg)");
    xlim([0,600]);
    subplot_ind = subplot_ind + 1;

    % % % controls % % %
    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind); 
    plot(show_t, ctrls(1,show_ind),"r","LineWidth",1.2); hold on;
    plot(show_t, ctrls(2,show_ind),"b","LineWidth",1.2); grid on;
    xlabel("Time since Release (s)"); ylabel("TVC");
    xlim([0,600]); ylim([-1,1]); hold off;
    subplot_ind = subplot_ind + 1;

    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind);
    plot(show_t, ctrls(3,show_ind)*100,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("Throttle (\%)");
    xlim([0,600]); ylim([0,100]);
    subplot_ind = subplot_ind + 1;

    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind);
    plot(show_t, ctrls(4,show_ind),"LineWidth",1.2); hold on; 
    plot(show_t, ctrls(5,show_ind),"LineWidth",1.2); grid on;
    plot(show_t, ctrls(6,show_ind),"LineWidth",1.2);
    plot(show_t, ctrls(7,show_ind),"LineWidth",1.2);
    plot(show_t, ctrls(8,show_ind),"LineWidth",1.2);
    plot(show_t, ctrls(9,show_ind),"LineWidth",1.2);
    plot(show_t, ctrls(10,show_ind),"LineWidth",1.2);
    plot(show_t, ctrls(11,show_ind),"LineWidth",1.2);
    xlabel("Time since Release (s)"); ylabel("ACS (\%)"); hold off
    xlim([0,600]);
    subplot_ind = subplot_ind + 1;

    axes(subplot_ind) = subplot(subplot_size(1),subplot_size(2),subplot_ind);
    plot(show_t, ctrls(12,show_ind)*180/pi,"LineWidth",1.2); hold on; 
    plot(show_t, ctrls(13,show_ind)*180/pi,"LineWidth",1.2); grid on;
    plot(show_t, ctrls(14,show_ind)*180/pi,"LineWidth",1.2);
    plot(show_t, ctrls(15,show_ind)*180/pi,"LineWidth",1.2);
    xlabel("Time since Release (s)"); ylabel("GF gimbal (deg)"); hold off
    xlim([0,600]);
    subplot_ind = subplot_ind + 1;
end