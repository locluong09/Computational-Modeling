function loc = rho_avg(p)
    rho_sc = 42;
    loc = rho_sc*(1 + 3e-5*(p-14.7))^-1;
end
