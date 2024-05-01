function param = cp_est(param)
% % James Barrowman Equation for Finding CP
param.brm.L_N = param.l.plf; % length of nosecone, m
param.brm.d_R = 0.609; % diameter of cylindrical body, m
param.brm.d_F = 0.6095; % diameter of nosecone, m
param.brm.d = param.brm.d_F;
param.brm.L_T = 0.0762; % length of transition (3in->m)
param.brm.X_P = param.brm.L_N; % distance from tip of nose to front of transition, m
param.brm.R = 0.609/2; % Rocket body diameter, m
param.brm.S = 6*0.3556; % Length of fin, m (6x14 in)
param.brm.C_R = 0.0254; % Fin thickness, m (1 in)
param.brm.C_T = param.brm.C_R;
param.brm.L_F = param.brm.S; % Length of fin mid-chord line, m
param.brm.X_R = 0; 
param.brm.X_B = param.d.enge1st; % Distance from nose tip to fin leading edge, m
param.brm.N = 4; % Number of fins
param.brm.C_NN = 2;
param.brm.X_N = 0.466*param.brm.L_N; % Ogive shaped nose cone
param.brm.C_NT = 2*((param.brm.d_R/param.brm.d)^2-(param.brm.d_F/param.brm.d)^2);
param.brm.X_T = param.brm.X_P+param.brm.L_T/3*(1+1/(1+param.brm.d_F/param.brm.d_R));
param.brm.C_NF = (1+param.brm.R/(param.brm.S+param.brm.R))*...
    (4*param.brm.N*(param.brm.S/param.brm.d)^2)/...
    (1+sqrt(1+(2*param.brm.L_F/(param.brm.C_R+param.brm.C_T))^2));
param.brm.X_F = param.brm.X_B + param.brm.X_R/3*(param.brm.C_R+2*param.brm.C_T)/(param.brm.C_R+param.brm.C_T) +...
    1/6*((param.brm.C_R+param.brm.C_T)-(param.brm.C_R*param.brm.C_T)/(param.brm.C_R+param.brm.C_T));
rcp_1st = (param.brm.C_NN*param.brm.X_N+param.brm.C_NT*param.brm.X_T+param.brm.C_NF*param.brm.X_F)/...
    (param.brm.C_NN+param.brm.C_NT+param.brm.C_NF);
rcp_2nd = (param.brm.C_NN*param.brm.X_N+param.brm.C_NT*param.brm.X_T)/(param.brm.C_NN+param.brm.C_NT);
param.rcp_1st = [-rcp_1st; 0; 0];
param.rcp_2nd = [-rcp_2nd; 0; 0];
end