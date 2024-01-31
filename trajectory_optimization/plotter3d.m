%% FUNCTION - 3D Trajectory Plotter
function plotter3d(x,N_1st,N_2nd,param)
    N = N_1st + N_2nd;
    unsc_x = x(param.init_ind_ptr(1)+(1:N))*param.scales.length;
    unsc_y = x(param.init_ind_ptr(2)+(1:N))*param.scales.length;
    unsc_z = x(param.init_ind_ptr(3)+(1:N))*param.scales.length;
    unsc_vx = x(param.init_ind_ptr(4)+(1:N))*param.scales.speed;
    unsc_vy = x(param.init_ind_ptr(5)+(1:N))*param.scales.speed;
    unsc_vz = x(param.init_ind_ptr(6)+(1:N))*param.scales.speed;
    
    figure; plot3(unsc_x,unsc_y,unsc_z); hold on;
    quiver3(unsc_x,unsc_y,unsc_z,unsc_vx,unsc_vy,unsc_vz);
    xl = xlim; yl = ylim; zl = zlim; 
    [X,Y,Z] = sphere(100); 
    surf(X*param.scales.length,Y*param.scales.length,Z*param.scales.length,"EdgeColor","none");
    axis equal; xlim(xl); ylim(yl); zlim(zl);
    xlabel("x"); ylabel("y"); zlabel("z");
end