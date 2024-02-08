%% FUNCTION - obtain real(unscaled) atmosphere by altitudeplf_dropalti_ref
% from NASA's interpolation on temp/pressure with altitude
function profile = atmo(h)
    h_trop = h(h<10999);
    h_lstrat = h(h>=10999&h<24999);
    h_ustrat = h(h>=24999);

    T_trop = 15.04-0.00649*h_trop;
    T_lstrat = -56.46*ones(size(h_lstrat));
    T_ustrat = -131.21+0.00299*h_ustrat;

    P_trop = 101.29*((T_trop+273.15)/288.08).^(5.256);
    P_lstrat = 22.65*exp(1.73-0.000157*h_lstrat);
    P_ustrat = 2.488*((T_ustrat+273.15)/216.6).^(-11.388);

    profile.T = [T_trop T_lstrat T_ustrat]+273.15;
    profile.P = [P_trop P_lstrat P_ustrat]*1000;
    profile.rho = profile.P./(286.9*profile.T);
    profile.h = h;
