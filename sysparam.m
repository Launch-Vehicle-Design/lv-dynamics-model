function [param,envir] = sysparam()

    % param.m0 = 2267.962;
    param.lPLF = 18*0.0254;
    param.dia = 24*0.0254;
    param.S = pi*param.dia^2/4;

    param.mPL = 40;
    
    % 1st stage
    param.m0 = 1816.7142;
    param.ms1 = 155.8821;
    param.mp1 = 1402.9389;
    param.Isp1 = 369.5;
    param.tbyw1 = 1.6;
    param.l1 = 3;
    
    % 2nd stage
    param.m02 = 257.8932;
    param.ms2 = 25.5472;
    param.mp2 = 187.346;
    param.Isp2 = 369.5;
    param.tbyw2 = 0.8;
    param.l2 = 3.5;

    param.L = param.l1+param.l2+param.lPLF;

    envir.g0 = 9.80665;
    envir.Rp = 6378100;
    
