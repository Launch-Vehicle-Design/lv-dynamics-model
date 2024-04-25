function param = extract_bottomup(vehicle_sizing, param)
scales_mass = 1; scales_length = 1;
if isfield(param,"scales")
    scales_mass = param.scales.mass;
    scales_length = param.scales.length;
end

stg1st = vehicle_sizing.first_stage;
stg2nd = vehicle_sizing.second_stage;
fn1st = fieldnames(stg1st);
fn2nd = fieldnames(stg2nd);

propname1st = 'solidprop';
oxname2nd = 'ox';
fuelname2nd = 'fuel';
plname2nd = 'payload';

ind_prop_1st = strcmp(fn1st,propname1st);
ind_ox_2nd = strcmp(fn2nd,oxname2nd);
ind_fuel_2nd = strcmp(fn2nd,fuelname2nd);
ind_pl_2nd = strcmp(fn2nd,plname2nd);
ind_nonprop_1st = ~ind_prop_1st;
ind_nonprop_2nd = ~(ind_ox_2nd|ind_fuel_2nd|ind_pl_2nd);

% dim variable arranged as following in SI units
% [length, mass, Izz, apx thickness, CG (rela2 fwd edge of component), CG (rela2 nose of rocket)]
ms1 = 0; cgs1 = 0;
ms2 = 0; cgs2 = 0;

% 1st stage sizing
prop1st = stg1st.(fn1st{ind_prop_1st}).dims;
index_nonprop_1st = find(ind_nonprop_1st);
for i = 1:length(index_nonprop_1st)
    index = index_nonprop_1st(i);
    dims = stg1st.(fn1st{index}).dims;
    weighted_cg = cgs1*ms1;
    ms1 = ms1 + dims(2);
    cgs1 = (weighted_cg + dims(2)*dims(6))/ms1;
end
mp1 = prop1st(2); m1 = mp1 + ms1;
% propellant and structure mass correction for residuals 1%
percent = 0.01;
mp1 = mp1/(1+percent); ms1 = ms1 + mp1*percent - 0.2*stg1st.interstage.dims(2);
cg1 = (cgs1*ms1 + mp1*prop1st(6))/m1;

% 2nd stage sizing
ox2nd = stg2nd.(fn2nd{ind_ox_2nd}).dims;
fuel2nd = stg2nd.(fn2nd{ind_fuel_2nd}).dims;
pl2nd = stg2nd.(fn2nd{ind_pl_2nd}).dims;
index_nonprop_2nd = find(ind_nonprop_2nd);
for i = 1:length(index_nonprop_2nd)
    index = index_nonprop_2nd(i);
    dims = stg2nd.(fn2nd{index}).dims;
    weighted_cg = cgs2*ms2;
    ms2 = ms2 + dims(2);
    cgs2 = (weighted_cg + dims(2)*dims(6))/ms2;
end
mp2 = ox2nd(2)+fuel2nd(2); m2 = mp2 + ms2;
% propellant and structure mass correction for residuals 0.7% for controls
percent = 0.007;
mp2 = mp2/(1+percent); ms2 = ms2 + mp2*percent + 0.2*stg1st.interstage.dims(2);
cg2 = (cgs2*ms2 + ox2nd(2)*ox2nd(6) + fuel2nd(2)*fuel2nd(6))/m2;

mpl = pl2nd(2); m0 = m1 + m2 + mpl;
cg = (m2*cg2 + m1*cg1 + mpl*pl2nd(6))/m0;

param.m0 = m0/scales_mass;
param.mpl = mpl/scales_mass;
param.cg = cg;

param.m1 = m1/scales_mass;
param.ms1 = ms1/scales_mass;
param.mp1 = mp1/scales_mass;
param.cgs1 = cgs1/scales_length;
param.cg1 = cg1/scales_length;

param.m02 = (m2+mpl)/scales_mass;
param.m2 = m2/scales_mass;
param.ms2 = ms2/scales_mass;
param.mp2 = mp2/scales_mass;
param.cgs2 = cgs2/scales_length;
param.cg2 = cg2/scales_length;

param.dims_prop1st = prop1st;
param.dims_ox2nd = ox2nd;
param.dims_fuel2nd = fuel2nd;
param.dims_pl2nd = pl2nd;
end

