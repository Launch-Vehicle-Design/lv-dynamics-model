function param = extract_bottomup(vehicle_sizing, param)
visualize = false;

scales_mass = 1; scales_length = 1;
scales_density = scales_mass/scales_length^3;
if isfield(param,"scales")
    scales_mass = param.scales.mass;
    scales_length = param.scales.length;
    scales_density = param.scales.density;
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

mfn = fieldnames(m);
for i = 1:numel(mfn)
    m.(mfn{i}) = m.(mfn{i})/scales_mass;
end

% % % component lengths % % %
length_dim = 1;
l.prop1st = stg1st.(prop1st).dims(length_dim);
l.case1st = stg1st.(case1st_solidcasing).dims(length_dim);
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
l.plpaf = stg2nd.(pl).dims(length_dim);

lfn = fieldnames(l);
for i = 1:numel(lfn)
    l.(lfn{i}) = l.(lfn{i})/scales_length;
end

% % % component top from the nose % % %
lcg_dim = 5; gcg_dim = 6;
d.prop1st = stg1st.(prop1st).dims(gcg_dim)-stg1st.(prop1st).dims(lcg_dim);
d.case1st = stg1st.(case1st_solidcasing).dims(gcg_dim)-stg1st.(case1st_solidcasing).dims(lcg_dim);
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
d.plpaf = stg2nd.(pl).dims(gcg_dim)-stg2nd.(pl).dims(lcg_dim);

dfn = fieldnames(d);
for i = 1:numel(dfn)
    d.(dfn{i}) = d.(dfn{i})/scales_length;
end

% % % dimension constants % % %
univsl_diam = 23.5*0.0254/scales_length;
enge1st_diam = 10.9*0.0254/scales_length;
enge2nd_diam = 14.9*0.0254/scales_length;
solid_density = 1783.8/scales_density;
ox_density = 1448/scales_density;
fuel_density = 1000/scales_density;
oxtank_length = 0.4376/scales_length;
fueltank_length = 0.5/scales_length;
pl_size = 1*0.3048/scales_length;

% % % component visualizer % % %
if visualize
figure;
location_shift = @(shape,d) shape + [d*ones([1,size(shape,2)]); zeros([1,size(shape,2)])];
% payload fairing
shape_plot = [0 -l.plf -l.plf 0; 0 univsl_diam/2 -univsl_diam/2 0];
shape_plot = location_shift(shape_plot,-d.plf);
plot(shape_plot(1,:),shape_plot(2,:),"k","LineWidth",1.2); hold on
% payload with payload attachment fitting
shape_plot = [0 0 -l.plpaf -l.plpaf 0 0; 0 pl_size/2 pl_size/2 -pl_size/2 -pl_size/2 0];
shape_plot = location_shift(shape_plot,-d.plpaf);
plot(shape_plot(1,:),shape_plot(2,:),"b","LineWidth",1.2); grid on
% second stage liquid tanks
shape_plot = [0 0 -fueltank_length-oxtank_length -fueltank_length-oxtank_length -fueltank_length -fueltank_length -fueltank_length 0 0;
    0 univsl_diam/2 univsl_diam/2 -univsl_diam/2 -univsl_diam/2 univsl_diam/2 -univsl_diam/2 -univsl_diam/2 0];
shape_plot = location_shift(shape_plot,-d.case2nd);
plot(shape_plot(1,:),shape_plot(2,:),"r","LineWidth",1.2); axis equal
% second stage engine
shape_plot = [0 -l.enge2nd -l.enge2nd 0; 0 enge2nd_diam/2 -enge2nd_diam/2 0];
shape_plot = location_shift(shape_plot,-d.enge2nd);
plot(shape_plot(1,:),shape_plot(2,:),"m","LineWidth",1.2);
% interstage
shape_plot = [0 0 -l.interstg -l.interstg 0 0; 0 univsl_diam/2 univsl_diam/2 -univsl_diam/2 -univsl_diam/2 0];
shape_plot = location_shift(shape_plot,-d.interstg);
plot(shape_plot(1,:),shape_plot(2,:),"k","LineWidth",1.2);
% first stage solid casing
shape_plot = [0 0 -l.case1st -l.case1st 0 0; 0 univsl_diam/2 univsl_diam/2 -univsl_diam/2 -univsl_diam/2 0];
shape_plot = location_shift(shape_plot,-d.case1st);
plot(shape_plot(1,:),shape_plot(2,:),"b","LineWidth",1.2);
% first stage aft skirt
shape_plot = [0 0 -l.afts1st -l.afts1st 0 0; 0 univsl_diam/2 univsl_diam/2 -univsl_diam/2 -univsl_diam/2 0];
shape_plot = location_shift(shape_plot,-d.afts1st);
plot(shape_plot(1,:),shape_plot(2,:),"r","LineWidth",1.2);
% first stage engine
shape_plot = [0 -l.enge1st -l.enge1st 0; 0 enge1st_diam/2 -enge1st_diam/2 0];
shape_plot = location_shift(shape_plot,-d.enge1st);
plot(shape_plot(1,:),shape_plot(2,:),"m","LineWidth",1.2);
end

% % % CG from the nose % % %
rcg.case1st = [-d.case1st; 0; 0] + funcs.rc_1ststg_pt(univsl_diam,l.case1st);
rcg.enge1st = [-d.enge1st; 0; 0] + funcs.rc_1ststg_eng(enge1st_diam,l.enge1st);
rcg.interstg = [-d.interstg; 0; 0] + funcs.rc_intstg(univsl_diam,l.interstg);
rcg.afts1st = [-d.afts1st; 0; 0] + funcs.rc_aft(univsl_diam,l.afts1st);
rcg.gnc1st = [-d.gnc1st; 0; 0] + funcs.rc_1stgnc(univsl_diam,l.gnc1st);

rcg.case2nd = [-d.case2nd; 0; 0] + funcs.rc_2ndstg_pt(univsl_diam,oxtank_length,fueltank_length);
rcg.enge2nd = [-d.enge2nd; 0; 0] + funcs.rc_2ndstg_eng(enge1st_diam,l.enge2nd);
rcg.gnc2nd = [-d.gnc2nd; 0; 0] + funcs.rc_2ndgnc(univsl_diam,l.gnc2nd);
rcg.plpaf = [-d.plpaf; 0; 0] + funcs.rc_pl(pl_size);
rcg.plf = [-d.plf; 0; 0] + funcs.rc_plf(univsl_diam,l.plf);

% % % fixed structure components MOI to nose cone % % %
% first step structure
m_1st = m.gnc1st + m.afts1st + m.interstg + m.case1st + m.enge1st;
cg_1st = (rcg.gnc1st*m.gnc1st + rcg.afts1st*m.afts1st + rcg.interstg*m.interstg + ...
    rcg.case1st*m.case1st + rcg.enge1st*m.enge1st)/m_1st;
moinc_1st = funcs.moic_1stgnc(m.gnc1st,univsl_diam,l.gnc1st) + funcs.moia(m.gnc1st,-rcg.gnc1st) + ...
    funcs.moic_aft(m.afts1st,univsl_diam,l.afts1st) + funcs.moia(m.afts1st,-rcg.afts1st) + ...
    funcs.moic_intstg(m.interstg,univsl_diam,l.interstg) + funcs.moia(m.interstg,-rcg.interstg) + ...
    funcs.moic_1ststg_pt(m.case1st,univsl_diam,l.case1st) + funcs.moia(m.case1st,-rcg.case1st) + ...
    funcs.moic_1ststg_eng(m.enge1st,enge1st_diam,l.enge1st) + funcs.moia(m.enge1st,-rcg.enge1st);
% second step structure
m_2nd = m.gnc2nd + m.case2nd + m.enge2nd;
cg_2nd = (rcg.gnc2nd*m.gnc2nd + rcg.case2nd*m.case2nd + rcg.enge2nd*m.enge2nd)/m_2nd;
moinc_2nd = funcs.moic_2ndgnc(m.gnc2nd,univsl_diam,l.gnc2nd) + funcs.moia(m.gnc2nd,-rcg.gnc2nd) + ...
    funcs.moic_2ndstg_pt(m.case2nd,univsl_diam,oxtank_length,fueltank_length) + funcs.moia(m.case2nd,-rcg.enge2nd) + ...
    funcs.moic_2ndstg_eng(m.enge2nd,enge2nd_diam,l.enge2nd) + funcs.moia(m.enge2nd,-rcg.gnc2nd);
% payload and payload attachment fitting
moinc_pl = funcs.moic_pl(m.plpaf,pl_size) + funcs.moia(m.plpaf,[d.plpaf; 0; 0]-funcs.rc_pl(pl_size));
% payload fairing
moinc_plf = funcs.moic_plf(m.plf,univsl_diam,l.plf) + funcs.moia(m.plf,[d.plf; 0; 0]-funcs.rc_plf(univsl_diam,l.plf));

% % % first stage fixed structure MOI wrt CGtotal % % %
m_1ststg = m_1st + m_2nd + m.plpaf + m.plf;
cg_1ststg = (m_1st*cg_1st + m_2nd*cg_2nd + m.plpaf*rcg.plpaf + m.plf*rcg.plf)/m_1ststg;
moinc_1ststg = moinc_1st + moinc_2nd + moinc_pl + moinc_plf;
moic_1ststg = moinc_1ststg - funcs.moia(m_1ststg,-cg_1ststg);

% % % second stage fixed structure MOI wrt CG2 with PLF % % %
m_2ndstg_wplf = m_2nd + m.plpaf + m.plf;
cg_2ndstg_wplf = (m_2nd*cg_2nd + m.plpaf*rcg.plpaf + m.plf*rcg.plf)/m_2ndstg_wplf;
moinc_2ndstg_wplf = moinc_2nd + moinc_pl + moinc_plf;
moic_2ndstg_wplf = moinc_2ndstg_wplf - funcs.moia(m_2ndstg_wplf,-cg_2ndstg_wplf);

% % % second stage fixed structure MOI wrt CG2 w/o PLF % % %
m_2ndstg = m_2nd + m.plpaf;
cg_2ndstg = (m_2nd*cg_2nd + m.plpaf*rcg.plpaf)/m_2ndstg;
moinc_2ndstg = moinc_2nd + moinc_pl;
moic_2ndstg = moinc_2ndstg - funcs.moia(m_2ndstg,-cg_2ndstg);

% % % % % pass parameters % % % % %

param.m0 = m0/scales_mass;
param.mpl = mpl/scales_mass;
param.cg = cg/scales_length;

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

param.funcs = funcs;

param.m_1ststg = m_1ststg;
param.m_2ndstg_wplf = m_2ndstg_wplf;
param.m_2ndstg = m_2ndstg;

param.moic_1ststg = moic_1ststg;
param.moic_2ndstg_wplf = moic_2ndstg_wplf;
param.moic_2ndstg = moic_2ndstg;

param.moinc_1ststg = moinc_1ststg;
param.moinc_2ndstg_wplf = moinc_2ndstg_wplf;
param.moinc_2ndstg = moinc_2ndstg;

param.cg_1ststg = cg_1ststg;
param.cg_2ndstg_wplf = cg_2ndstg_wplf;
param.cg_2ndstg = cg_2ndstg;

param.m = m;
param.l = l;
param.d = d;
param.rcg = rcg;
param.univsl_diam = univsl_diam;
param.solid_density = solid_density;
param.ox_density = ox_density;
param.fuel_density = fuel_density;

%% DEFINE PROPELLANT MOI FUNCTIONS
% % % first stage consuming propellant - end burner geometry % % %
param.rnc_prop1st = @(m) [-d.prop1st; 0; 0]+funcs.rc_1ststg_p(param.univsl_diam,m/param.solid_density/(pi*(param.univsl_diam/2)^2));
param.drnc_prop1st = @(m,dm) funcs.drc_1ststg_p(param.univsl_diam,m/param.solid_density/(pi*(param.univsl_diam/2)^2),0,dm/param.solid_density/(pi*(param.univsl_diam/2)^2));
param.moinc_prop1st = @(m) param.funcs.moic_1ststg_p_rgb(m,param.univsl_diam,m/param.solid_density/(pi*(param.univsl_diam/2)^2)) + ...
    param.funcs.moia(m,param.rnc_prop1st(m));
param.dmoinc_prop1st = @(m,dm) param.funcs.dmoic_1ststg_p_rbg(m,param.univsl_diam,m/param.solid_density/(pi*(param.univsl_diam/2)^2),dm,0,dm/param.solid_density/(pi*(param.univsl_diam/2)^2)) + ...
    param.funcs.dmoia(m,param.rnc_prop1st(m),dm,param.drnc_prop1st(m,dm));

% % % second stage consuming propellant - pendulum model % % %
param.rnc_ox2nd = @(m) [-d.ox2nd; 0; 0]+funcs.rc_2ndstg_p(param.univsl_diam,m/param.ox_density/(pi*(param.univsl_diam/2)^2));
param.rnc_fuel2nd = @(m) [-d.fuel2nd; 0; 0]+funcs.rc_2ndstg_p(param.univsl_diam,m/param.fuel_density/(pi*(param.univsl_diam/2)^2));
param.drnc_ox2nd = @(m,dm) funcs.drc_2ndstg_p(param.univsl_diam,m/param.ox_density/(pi*(param.univsl_diam/2)^2),0,dm/param.ox_density/(pi*(param.univsl_diam/2)^2));
param.drnc_fuel2nd = @(m,dm) funcs.drc_2ndstg_p(param.univsl_diam,m/param.fuel_density/(pi*(param.univsl_diam/2)^2),0,dm/param.fuel_density/(pi*(param.univsl_diam/2)^2));
param.moinc_ox2nd = @(m) param.funcs.I0(m,param.univsl_diam,m/param.ox_density/(pi*(param.univsl_diam/2)^2)) + param.funcs.moia(m,param.rnc_ox2nd(m));
param.moinc_fuel2nd = @(m) param.funcs.I0(m,param.univsl_diam,m/param.fuel_density/(pi*(param.univsl_diam/2)^2)) + ...
    param.funcs.moia(m,param.rnc_fuel2nd(m));
param.dmoinc_ox2nd = @(m,dm) param.funcs.dI0(m,param.univsl_diam,m/param.ox_density/(pi*(param.univsl_diam/2)^2),dm,0,dm/param.ox_density/(pi*(param.univsl_diam/2)^2)) + ...
    param.funcs.dmoia(m,param.rnc_ox2nd(m),dm,param.drnc_ox2nd(m,dm));
param.dmoinc_fuel2nd = @(m,dm) param.funcs.dI0(m,param.univsl_diam,m/param.fuel_density/(pi*(param.univsl_diam/2)^2),dm,0,dm/param.fuel_density/(pi*(param.univsl_diam/2)^2)) + ...
    param.funcs.dmoia(m,param.rnc_fuel2nd(m),dm,param.drnc_fuel2nd(m,dm));

% % % total vehicle MOI first stage % % %
param.rc1st = @(m,mox,mfuel) (param.cg_1ststg*param.m_1ststg + ...
    param.rnc_prop1st(param.mp1-(param.m0-m-(param.m.ox2nd-mox)-(param.m.fuel2nd-mfuel)))*(param.mp1-(param.m0-m-(param.m.ox2nd-mox)-(param.m.fuel2nd-mfuel))) + ...
    param.rnc_ox2nd(mox)*mox + param.rnc_fuel2nd(mfuel)*mfuel)/m;
param.moic1st = @(m,mox,mfuel) param.moinc_1ststg + ...
    param.moinc_prop1st(param.mp1-(param.m0-m-(param.m.ox2nd-mox)-(param.m.fuel2nd-mfuel))) + ...
    param.moinc_ox2nd(mox) + param.moinc_fuel2nd(mfuel) - ...
    param.funcs.moia(m,param.rc1st(m,mox,mfuel));
param.dmoic1st = @(m,mox,mfuel,dm,dmox,dmfuel) param.dmoinc_prop1st(param.mp1-(param.m0-m-(param.m.ox2nd-mox)-(param.m.fuel2nd-mfuel)),dm-dmox-dmfuel) + ...
    param.dmoinc_ox2nd(mox,dmox) + param.dmoinc_fuel2nd(mfuel,dmfuel) - ...
    param.funcs.dmoia(m,param.rc1st(m,mox,mfuel),dm,0); % ignore the center of gravity shift

% % % total vehicle MOI second stage with PLF % % %
param.rc2nd_wplf = @(m,mox,mfuel) (param.cg_2ndstg_wplf*param.m_2ndstg_wplf + ...
    param.rnc_ox2nd(mox)*mox + param.rnc_fuel2nd(mfuel)*mfuel)/m;
param.moic2nd_wplf = @(m,mox,mfuel) param.moinc_2ndstg_wplf + ...
    param.moinc_ox2nd(mox) + param.moinc_fuel2nd(mfuel) - ...
    param.funcs.moia(m,param.rc2nd_wplf(m,mox,mfuel));
param.dmoic2nd_wplf = @(m,mox,mfuel,dm,dmox,dmfuel) param.dmoinc_ox2nd(mox,dmox) + ...
    param.dmoinc_fuel2nd(mfuel,dmfuel) - ...
    param.funcs.dmoia(m,param.rc2nd_wplf(m,mox,mfuel),dm,0); % ignore the center of gravity shift

% % % total vehicle MOI second stage w/o PLF % % %
param.rc2nd = @(m,mox,mfuel) (param.cg_2ndstg*param.m_2ndstg + ...
    param.rnc_ox2nd(mox)*mox + param.rnc_fuel2nd(mfuel)*mfuel)/m;
param.moic2nd = @(m,mox,mfuel) param.moinc_2ndstg + ...
    param.moinc_ox2nd(mox) + param.moinc_fuel2nd(mfuel) - ...
    param.funcs.moia(m,param.rc2nd(m,mox,mfuel));
param.dmoic2nd = param.dmoic2nd_wplf; % ignore the center of gravity shift
end

