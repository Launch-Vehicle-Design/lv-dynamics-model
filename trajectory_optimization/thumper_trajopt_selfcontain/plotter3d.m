%% FUNCTION - 3D Trajectory Plotter
function plotter3d(x,param,new_fig)
    if exist("new_fig","var") && new_fig
        figure; end
    N_1st = param.N_1st; N_2nd = param.N_2nd;
    N = N_1st + N_2nd;
    unsc_x = x(param.init_ind_ptr(1)+(1:N))*param.scales.length;
    unsc_y = x(param.init_ind_ptr(2)+(1:N))*param.scales.length;
    unsc_z = x(param.init_ind_ptr(3)+(1:N))*param.scales.length;
    unsc_vx = x(param.init_ind_ptr(4)+(1:N))*param.scales.speed;
    unsc_vy = x(param.init_ind_ptr(5)+(1:N))*param.scales.speed;
    unsc_vz = x(param.init_ind_ptr(6)+(1:N))*param.scales.speed;
    
    show_num = 30; show_ind = 1:N/show_num:N;
    plot3(unsc_x,unsc_y,unsc_z,"LineWidth",2,"Color","r"); hold on;
    quiver3(unsc_x(show_ind),unsc_y(show_ind),unsc_z(show_ind),unsc_vx(show_ind),unsc_vy(show_ind),unsc_vz(show_ind),"LineWidth",3,"Color","b");
    xl = xlim; yl = ylim; zl = zlim; 
    [X,Y,Z] = sphere(100); 
    surf(X*param.scales.length,Y*param.scales.length,Z*param.scales.length,"EdgeColor","none","FaceColor",[135 206 235]/255,"FaceAlpha",0.5,"FaceLighting","gouraud");
    axis equal; xlim(xl); ylim(yl); zlim(zl);
    xlabel("x"); ylabel("y"); zlabel("z");
end