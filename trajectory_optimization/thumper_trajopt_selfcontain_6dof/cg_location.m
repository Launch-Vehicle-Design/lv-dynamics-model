
load("thumper_straj_cg3dof.mat");
param = sysparam(log_param,log_param.scales);

mass = states(7,:);
rcg = param.rc1st(1,param.m.ox2nd,param.m.fuel2nd)*ones(size(mass));
% stage 1 center of mass
for i = 1:length(stg1_ind)
    ind = stg1_ind(i);
    curr_mass = mass(ind)/param.scales.mass;
    rcg(:,ind) = param.rc1st(curr_mass,param.m.ox2nd,param.m.fuel2nd);
end
% stage separation center of mass
rcg(:,stg1p5_ind) = param.rc2nd(param.m02,param.m.ox2nd,param.m.fuel2nd)*ones(1,length(stg1p5_ind));

% stage 2 center of mass
for i = 1:length(stg2_ind)
    ind = stg2_ind(i);
    curr_mass = mass(ind)/param.scales.mass;
    prop_mass_burned = param.m02 - curr_mass;
    fuel_mass = prop_mass_burned/(param.OFratio+1);
    ox_mass = fuel_mass*param.OFratio;
    rcg(:,ind) = param.rc2nd(curr_mass,param.m.ox2nd-ox_mass,param.m.fuel2nd-fuel_mass);
end

rcg = rcg*param.scales.length;