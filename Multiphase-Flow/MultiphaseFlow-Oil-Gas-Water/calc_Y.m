function out = calc_X(k1, k2, h1, h2, delta_x1, delta_x2, delta_y1, delta_y2)
    beta_c = 1.127e-3;
    %geometric factor
    Ax1 = delta_y1*h1;
    Ax2 = delta_y2*h2;
    out = 2*beta_c/(delta_x1/(Ax1*k1) + delta_x2/(Ax2*k2));
end