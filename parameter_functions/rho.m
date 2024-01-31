function r = rho(h)
    h0 = 7640;
    rho0 = 1.225;
    r = rho0.*exp(-h/h0);