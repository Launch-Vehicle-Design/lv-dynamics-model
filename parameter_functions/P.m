function p = P(h)
    h0 = 7640;
    p0 = 101325;
    p = p0.*exp(-h/h0);