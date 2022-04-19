function y = tran_Y(k1, k2, h1, h2, miu, B, dx1, dx2, dy1, dy2 )
    beta_c = 1.127e-3;
    Ay1 = dx1*h1;
    Ay2 = dx2*h2;
    y = 2*beta_c/(miu*B)/(dy1/(Ay1*k1)+dy2/(Ay2*k2));
end