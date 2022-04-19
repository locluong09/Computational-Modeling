function out = calcGamma(delta_x, delta_y, h)
    delta_t = 5;
    phi = 0.2;
    c = 3e-6;
    alpha_c = 5.615;
    out = delta_x*delta_y*h*phi/(alpha_c*delta_t);
end
