function out = Jacobian(P_oil, Sg, Sw, Order, Pre_init, Sg_init, Sw_init)
    [Row, Col] = size(Order);
    m = Row - 2;
    n = Col - 2;
    x = zeros(m*n*3,1);
    for i = 1:Row
        for j = 1:Col
            if Order(i,j) > 0
                x(Order(i,j)*3-2) = P_oil(i,j);
                x(Order(i,j)*3-1) = Sg(i,j);
                x(Order(i,j)*3) = Sw(i,j);
            end
        end
    end
    dp = 1e-3;
    dSg = 1e-7;
    dSw = 1e-7;

    J = zeros(27,27);
    [T,Q] = construct_Residual_matrix(P_oil, Sg, Sw, Order, Pre_init, Sg_init, Sw_init);
    Residual = T*x - Q;
    for i = 1:Row
        for j = 1:Col
            if Order(i,j) > 0
                P_oil_new = P_oil;
                P_oil_new(i,j) = P_oil_new(i,j) +  dp;
                

                x_new1 = x;
                x_new1(Order(i,j)*3-2) = x_new1(Order(i,j)*3-2) + dp;

                [To,Qo] = construct_Residual_matrix(P_oil_new, Sg, Sw, Order, Pre_init, Sg_init, Sw_init);
                R_new1 = To*x_new1 - Qo;
                J(:,Order(i,j)*3-2) = (R_new1 - Residual)/dp;
               
                Sg_new = Sg;
                Sg_new(i,j) = Sg_new(i,j) + dSg;

                x_new2 = x;
                x_new2(Order(i,j)*3-1) = x_new2(Order(i,j)*3-1)+dSg;

                [Tg,Qg] = construct_Residual_matrix(P_oil, Sg_new, Sw, Order, Pre_init, Sg_init, Sw_init);
                R_new2 = Tg*x_new2 - Qg;
                
                J(:,Order(i,j)*3-1) = (R_new2 - Residual)/dSg;


                Sw_new = Sw;
                Sw_new(i,j) = Sw_new(i,j) + dSw;

                x_new3 = x;
                x_new3(Order(i,j)*3) = x_new3(Order(i,j)*3) + dSw;

                [Tw,Qw] = construct_Residual_matrix(P_oil, Sg, Sw_new, Order, Pre_init, Sg_init, Sw_init);
                R_new3 = Tw*x_new3 - Qw;
               
                J(:,Order(i,j)*3) = (R_new3 - Residual)/dSw;
            end


        end
    end

    out = J;
end