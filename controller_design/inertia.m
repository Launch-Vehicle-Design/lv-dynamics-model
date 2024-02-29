function iner = inertia(m,dm,param)
pl2nd = param.dims_pl2nd;
if m > param.m0 - param.mp1
    mprop = m-(param.m0-param.mp1);
    prop1st = param.dims_prop1st;
    cgp1 = mprop/param.mp1*prop1st(5) + prop1st(6)-prop1st(5);
    iner.cg = (param.mpl*pl2nd(6) + param.cg2*param.m2 + param.cgs1*param.ms1 + mprop*cgp1)/m;
elseif m > param.m0 - param.m1
    iner.cg = (param.mpl*pl2nd(6) + param.cg2*param.m2 + param.cgs1*param.ms1)/(param.mpl+param.m2+param.ms1);
elseif m > param.m02 - param.mp2
    mprop = m-(param.m02-param.mp2);
    ox2nd = param.dims_ox2nd; fuel2nd = param.dims_fuel2nd;
    mfuel = mprop/(1+param.OFratio); mox = mfuel*param.OFratio;
    cgf2 = mfuel/fuel2nd(2)*fuel2nd(5) + fuel2nd(6)-fuel2nd(5);
    cgo2 = mox/ox2nd(2)*ox2nd(5) + ox2nd(6)-ox2nd(5);
    iner.cg = (param.mpl*pl2nd(6) + param.cgs2*param.ms2 + mfuel*cgf2 + mox*cgo2)/m;
else
    iner.cg = (param.mpl*pl2nd(6) + param.cgs2*param.ms2)/(param.mpl+param.ms2);
end
end
