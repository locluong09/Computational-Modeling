function y = trans_X(k1, k2, h1, h2, miu, B, dx1, dx2, dy1, dy2)
    beta_c = 1.127e-3;
    Ax1 = dy1*h1;
    Ax2 = dy2*h2;
    y = 2*beta_c/(miu*B)/(dx1/(Ax1*k1)+dx2/(Ax2*k2));
end
