function axes = plotter3d(ind,state,traj,param,mpc_mode,mpc_states)
% % 3D visualization of launch vehicle from LVLH frame, centered at COM
ref_hrzn = 10;
if ~exist("mpc_mode","var") || ~exist("mpc_states","var") 
    mpc_mode = false;
end

% construct the LVLH frame
r = state(1:3); k = [0;0;1];
e_lvz = r/vecnorm(r);
e_esx = cross(k,e_lvz)/vecnorm(cross(k,e_lvz));
e_nry = cross(e_lvz,e_esx);
% rotation matrix from tcvh to iner
C_enu = [e_esx e_nry e_lvz];

% extract body frame
C_body = EP2C(state(7:10));
C_b2enu = C_enu'*C_body;

% vehicle graph plot
m = state(14); mox = state(15); mfuel = state(16);
lvgraph(lvrot(C_b2enu,param.body,param.rc1st(m,mox,mfuel))); hold on

% extract reference trajectory
low_ind = max(1,ind-ref_hrzn/2);
upp_ind = min(ind+ref_hrzn/2,size(traj.t,2));
ref_traj = traj.states(:,low_ind:upp_ind); 
ref_pos_enu = C_enu'*(ref_traj(1:3,:)-state(1:3));

% visualize the body frame
scale = 5; C_scaled = scale*C_b2enu;
quiver3(0,0,0,C_scaled(1,1),C_scaled(2,1),C_scaled(3,1),"Color",[1 0 0],"LineWidth",2);
quiver3(0,0,0,C_scaled(1,2),C_scaled(2,2),C_scaled(3,2),"Color",[0 1 0],"LineWidth",2);
quiver3(0,0,0,C_scaled(1,3),C_scaled(2,3),C_scaled(3,3),"Color",[0 0 1],"LineWidth",2);

% visualize control prediction horizon
if mpc_mode
    mpc_pos_enu = C_enu'*(mpc_states(1:3,:)-state(1:3));
    plot3(mpc_pos_enu(1,:),mpc_pos_enu(2,:),mpc_pos_enu(3,:),"b--","LineWidth",2);
end

% plot the reference trajectory at some horizon
plot3(ref_pos_enu(1,:),ref_pos_enu(2,:),ref_pos_enu(3,:),"Color",[0.7 0.7 0.7],"LineWidth",2);
axis equal; ylim([-100, 100]); hold off
end
