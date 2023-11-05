function gh = g(h)
    [~,envir] = sysparam();
    g0 = envir.g0;
    Rp = envir.Rp;
    gh = g0/(1+h/Rp)^2;