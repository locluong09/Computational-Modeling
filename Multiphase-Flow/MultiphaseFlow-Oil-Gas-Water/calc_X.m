function out = calc_Y(k1, k2, h1, h2, delta_x1, delta_x2, delta_y1, delta_y2)
    beta_c = 1.127e-3;
    Ay1 = delta_x1*h1;
    Ay2 = delta_x2*h2;
    out = 2*beta_c/(delta_y1/(Ay1*k1) + delta_y2/(Ay2*k2));
end