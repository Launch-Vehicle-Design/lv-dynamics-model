%% FUNCTION - State and Control Plotter
function axes = plotter(x,N_1st,N_2nd,param)
    N = N_1st + N_2nd;
    ntot = param.nstate + param.nctrl;
    t_1st = 0:x(ntot*N+1):x(ntot*N+1)*(N_1st-1);
    t_2nd = x(ntot*N+1)*N_1st:x(ntot*N+2):x(ntot*N+1)*N_1st+x(ntot*N+2)*(N_2nd-1);
    if isfield(param,"axes") axes = param.axes; end
    unsc_t = [t_1st t_2nd]*param.scales.time;
    unsc_x = x(param.init_ind_ptr(1)+(1:N))*param.scales.length;
    unsc_y = x(param.init_ind_ptr(2)+(1:N))*param.scales.length;
    unsc_z = x(param.init_ind_ptr(3)+(1:N))*param.scales.length;
    unsc_vx = x(param.init_ind_ptr(4)+(1:N))*param.scales.speed;
    unsc_vy = x(param.init_ind_ptr(5)+(1:N))*param.scales.speed;
    unsc_vz = x(param.init_ind_ptr(6)+(1:N))*param.scales.speed;
    unsc_ux = x(param.init_ind_ptr(8)+(1:N));
    unsc_uy = x(param.init_ind_ptr(9)+(1:N));
    rn = vecnorm([unsc_x';unsc_y';unsc_z'])'; vn = vecnorm([unsc_vx';unsc_vy';unsc_vz'])';
    pitch_ang = (unsc_vx.*unsc_x + unsc_vy.*unsc_y + unsc_vz.*unsc_z)./rn./vn;
    
    axes(1) = subplot(2,3,1);
    plot(unsc_t,vecnorm([unsc_x unsc_y unsc_z]')-param.earthR*param.scales.length,"k","LineWidth",1.2); grid on;
    % plot(unsc_t, unsc_x,"r","LineWidth",1.2); grid on; hold on;
    % plot(unsc_t, unsc_y,"b","LineWidth",1.2); plot(unsc_t, unsc_z,"k","LineWidth",1.2);
    xlabel("Time since Release (s)"); ylabel("Vehicle Position (m)"); hold off;

    axes(2) = subplot(2,3,2);
    plot(unsc_t,vecnorm([unsc_vx unsc_vy unsc_vz]'),"k","LineWidth",1.2); grid on;
    % plot(unsc_t, unsc_vx,"r","LineWidth",1.2); grid on; hold on;
    % plot(unsc_t, unsc_vy,"b","LineWidth",1.2); plot(unsc_t, unsc_vz,"k","LineWidth",1.2);
    xlabel("Time since Release (s)"); ylabel("Vehicle Velocity (m/s)"); hold off;

    axes(3) = subplot(2,3,3); 
    plot(unsc_t, x(param.init_ind_ptr(7)+(1:N))*param.scales.mass,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("Vehicle Mass (kg)");

    axes(4) = subplot(2,3,4); 
    plot(unsc_t, unsc_ux*180/pi,"r","LineWidth",1.2); grid on; hold on;
    plot(unsc_t, unsc_uy*180/pi,"b","LineWidth",1.2);
    xlabel("Time since Release (s)"); ylabel("Gimbal Angle (deg)"); hold off;

    axes(5) = subplot(2,3,5);
    plot(unsc_t,(pi/2-acos(pitch_ang))*180/pi,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("pitch, undirectional (deg)");

    axes(6) = subplot(2,3,6); 
    plot(unsc_t, x(param.init_ind_ptr(10)+(1:N))*100,"k","LineWidth",1.2); grid on
    xlabel("Time since Release (s)"); ylabel("Throttle (\%)");
end