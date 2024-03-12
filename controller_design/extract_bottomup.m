function param = extract_bottomup(vehicle_sizing, param)
scales_mass = 1; scales_length = 1;
if isfield(param,"scales")
    scales_mass = param.scales.mass;
    scales_length = param.scales.length;
end
funcs = moi_func;

stg1st = vehicle_sizing.first_stage;
stg2nd = vehicle_sizing.second_stage;
fn1st = fieldnames(stg1st);
fn2nd = fieldnames(stg2nd);

%% mass extraction
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
cg2 = (cgs2*ms2 + ox2nd(2)*ox2nd(6) + fuel2nd(2)*fuel2nd(6))/m2;

mpl = pl2nd(2); m0 = m1 + m2 + mpl;
cg = (m2*cg2 + m1*cg1 + mpl*pl2nd(6))/m0;

%% hard coded component list
prop1st = "solidprop";
case1st_insl = "insl"; 
case1st_solidcasing = "solidcasing";
enge1st = "engine";
afts1st = "aftskirt";
gnc1st = "avgnc";
interstg = "interstage";
plf = "PLF";

ox2nd = "ox";
fuel2nd = "fuel";
case2nd = "liqcasing";
enge2nd = "engine";
gnc2nd = "wiringgnc";
paf = "PAF";
pl = "payload";

% % % component masses % % %
mass_dim = 2;
m.prop1st = stg1st.(prop1st).dims(mass_dim);
m.case1st = stg1st.(case1st_insl).dims(mass_dim) + ...
    stg1st.(case1st_solidcasing).dims(mass_dim);
m.enge1st = stg1st.(enge1st).dims(mass_dim);
m.afts1st = stg1st.(afts1st).dims(mass_dim);
m.gnc1st = stg1st.(gnc1st).dims(mass_dim);
m.interstg = stg1st.(interstg).dims(mass_dim);
m.plf = stg1st.(plf).dims(mass_dim);

m.ox2nd = stg2nd.(ox2nd).dims(mass_dim);
m.fuel2nd = stg2nd.(fuel2nd).dims(mass_dim);
m.case2nd = stg2nd.(case2nd).dims(mass_dim);
m.enge2nd = stg2nd.(enge2nd).dims(mass_dim);
m.gnc2nd = stg2nd.(gnc2nd).dims(mass_dim);
m.plpaf = stg2nd.(pl).dims(mass_dim) + ...
    stg2nd.(paf).dims(mass_dim);

% % % component lengths % % %
length_dim = 1;
l.prop1st = stg1st.(prop1st).dims(length_dim);
l.case1st = stg1st.(case1st_insl).dims(length_dim) + ...
    stg1st.(case1st_solidcasing).dims(length_dim);
l.enge1st = stg1st.(enge1st).dims(length_dim);
l.afts1st = stg1st.(afts1st).dims(length_dim);
l.gnc1st = stg1st.(gnc1st).dims(length_dim);
l.interstg = stg1st.(interstg).dims(length_dim);
l.plf = stg1st.(plf).dims(length_dim);

l.ox2nd = stg2nd.(ox2nd).dims(length_dim);
l.fuel2nd = stg2nd.(fuel2nd).dims(length_dim);
l.case2nd = stg2nd.(gnc2nd).dims(length_dim);
l.enge2nd = stg2nd.(enge2nd).dims(length_dim);
l.gnc2nd = stg2nd.(gnc2nd).dims(length_dim);
l.plpaf = stg2nd.(pl).dims(length_dim) + ...
    stg2nd.(paf).dims(length_dim);

% % % component top from the nose % % %
lcg_dim = 5; gcg_dim = 6;
d.prop1st = stg1st.(prop1st).dims(gcg_dim)-stg1st.(prop1st).dims(lcg_dim);
d.case1st = stg1st.(case1st_insl).dims(gcg_dim)-stg1st.(case1st_insl).dims(lcg_dim) + ...
    stg1st.(case1st_solidcasing).dims(gcg_dim)-stg1st.(case1st_solidcasing).dims(lcg_dim);
d.enge1st = stg1st.(enge1st).dims(gcg_dim)-stg1st.(enge1st).dims(lcg_dim);
d.afts1st = stg1st.(afts1st).dims(gcg_dim)-stg1st.(afts1st).dims(lcg_dim);
d.gnc1st = stg1st.(gnc1st).dims(gcg_dim)-stg1st.(gnc1st).dims(lcg_dim);
d.interstg = stg1st.(interstg).dims(gcg_dim)-stg1st.(interstg).dims(lcg_dim);
d.plf = stg1st.(plf).dims(gcg_dim)-stg1st.(plf).dims(lcg_dim);

d.ox2nd = stg2nd.(ox2nd).dims(gcg_dim)-stg2nd.(ox2nd).dims(lcg_dim);
d.fuel2nd = stg2nd.(fuel2nd).dims(gcg_dim)-stg2nd.(fuel2nd).dims(lcg_dim);
d.case2nd = stg2nd.(gnc2nd).dims(gcg_dim)-stg2nd.(gnc2nd).dims(lcg_dim);
d.enge2nd = stg2nd.(enge2nd).dims(gcg_dim)-stg2nd.(enge2nd).dims(lcg_dim);
d.gnc2nd = stg2nd.(gnc2nd).dims(gcg_dim)-stg2nd.(gnc2nd).dims(lcg_dim);
d.plpaf = stg2nd.(pl).dims(gcg_dim)-stg2nd.(pl).dims(lcg_dim) + ...
    stg2nd.(paf).dims(gcg_dim)-stg2nd.(paf).dims(lcg_dim);

% % % moment of inertia % % %
univsl_diam = 23.5*0.0254;
enge1st_diam = 10.9*0.0254;
enge2nd_diam = 14.9*0.0254;
oxtank_length = 0.4376;
fueltank_length = 0.5;
pl_size = 1*0.3048;
% % % fixed structure components MOI to nose cone % % %
% first step structure
m_1st = m.gnc1st + m.afts1st + m.interstg + m.case1st + m.enge1st;
moinc_1st = funcs.moic_1stgnc(m.gnc1st,univsl_diam,l.gnc1st) + funcs.moia(m.gnc1st,[d.gnc1st; 0; 0]-funcs.rc_1stgnc(univsl_diam,l.gnc1st)) + ...
    funcs.moic_aft(m.afts1st,univsl_diam,l.afts1st) + funcs.moia(m.afts1st,[d.afts1st; 0; 0]-funcs.rc_aft(univsl_diam,l.afts1st)) + ...
    funcs.moic_intstg(m.interstg,univsl_diam,l.interstg) + funcs.moia(m.interstg,[d.interstg; 0; 0]-funcs.rc_intstg(univsl_diam,l.interstg)) + ...
    funcs.moic_1ststg_pt(m.case1st,univsl_diam,l.case1st) + funcs.moia(m.case1st,[d.case1st; 0; 0]-funcs.rc_1ststg_pt(univsl_diam,l.case1st)) + ...
    funcs.moic_1ststg_eng(m.enge1st,enge1st_diam,l.enge1st) + funcs.moia(m.enge1st,[d.enge1st; 0; 0]-funcs.rc_1ststg_eng(enge1st_diam,l.enge1st));
% second step structure
m_2nd = m.gnc2nd + m.case2nd + m.enge2nd;
moinc_2nd = funcs.moic_2ndgnc(m.gnc2nd,univsl_diam,l.gnc2nd) + funcs.moia(m.gnc2nd,[d.gnc2nd; 0; 0]-funcs.rc_2ndgnc(univsl_diam,l.gnc2nd)) + ...
    funcs.moic_2ndstg_pt(m.case2nd,univsl_diam,oxtank_length,fueltank_length) + funcs.moia(m.case2nd,[d.case2nd; 0; 0]-funcs.rc_2ndstg_pt) + ...
    funcs.moic_2ndstg_eng(m.enge2nd,enge2nd_diam,l.enge2nd) + funcs.moia(m.enge2nd,[d.enge2nd; 0; 0]-funcs.rc_2ndstg_eng(enge2nd_diam,l.enge2nd));
% payload and payload attachment fitting
moinc_pl = funcs.moic_pl(m.plpaf,pl_size) + funcs.moia(m.plpaf,[d.plpaf; 0; 0]-funcs.rc_pl(pl_size));
% payload fairing
moinc_plf = funcs.moic_plf(m.plf,univsl_diam,l.plf) + funcs.moia(m.plf,[d.plf; 0; 0]-funcs.rc_plf(univsl_diam,l.plf));

% % % first stage fixed structure MOI wrt CGtotal % % %
m_1ststg = m_1st + m_2nd + m.plpaf + m.plf;
moinc_1st = moinc_1st + moinc_2nd + moinc_pl + moinc_plf;

% % % second stage fixed structure MOI wrt CG2 with PLF % % %
m_2ndstg_wplf = m_2nd + m.plpaf + m.plf;
moinc_2nd_wplf = moinc_2nd + moinc_pl + moinc_plf;

% % % second stage fixed structure MOI wrt CG2 w/o PLF % % %
m_2ndstg = m_2nd + m.plpaf;
moinc_2nd = moinc_2nd + moinc_pl;

param.funcs = funcs;
param.nc2cg = @(moinc,m,d_cg) moinc - param.funcs.moia(m,[d_cg; 0; 0]);

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

