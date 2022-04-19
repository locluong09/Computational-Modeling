function y = B_avg(p)
    c_o = 3e-5;
    p_sc = 14.7;
    y = (1+c_o*(p-p_sc))^-1;
end
