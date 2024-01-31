function gh = g(h)
    g0 = 9.80665;
    Rp = 6378100;
    gh = g0./(1+h/Rp).^2;