function [out1, out2] = construct_Residual_matrix(P_oil, Sg, Sw, Order, P_init, Sg_init, Sw_init)
    [Row, Col] = size(Order);
    m = Row - 2;
    n = Col - 2;
    T = zeros((Row-2)*(Col-2)*3, (Row-2)*(Col-2)*3);
    Q = zeros(m*n*3,1);
    oil_prop  = load("oil_properties.txt");
    water_prop  = load("water_properties.txt");
    gas_prop  = load("gas_properties.txt");
    oil_water_rel  = load("oil_water_rel_table.txt");
    gas_oil_rel  = load("gas_oil_rel_table.txt");
    Perm = ones(3,3)*100;
    Thick = ones(3,3)*50;
    DX = ones(3,3)*500;
    DY = ones(3,3)*500;
    Phi = ones(3,3)*0.2;

    Perm = boundaries(Perm);
    Thickness = boundaries(Thick);
    DX = boundaries(DX);
    DY = boundaries(DY);
    poro = boundaries(Phi);
    por_new = zeros(Row,Col);
    %[Row, Col] = size(Perm);
    for i = 1:Row
        for j = 1:Col
            if Order(i,j) > 0
                Gamma = calcGamma(500,500,50);
                S_w = Sw(i,j);
                S_g = Sg(i,j);

                P_cow = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,4), S_w);
                P_cgo = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,4), S_g);

                B_o_rhs = linearInterpolation(oil_prop(:,1), oil_prop(:,3), P_oil(i,j));
                B_g_rhs = linearInterpolation(gas_prop(:,1), gas_prop(:,3), P_oil(i,j) + P_cgo);
                B_w_rhs = linearInterpolation(water_prop(:,1), water_prop(:,3), P_oil(i,j) - P_cow);

                R_so_rhs = linearInterpolation(oil_prop(:,1), oil_prop(:,5), P_oil(i,j));

                if Order(i-1,j) > 0
                    P_o = (P_oil(i-1,j) + P_oil(i,j))/2;
                    S_w_North = Sw(i-1,j);
                    S_g_North = Sg(i-1,j);
                    P_cow_North = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,4), S_w_North);
                    P_cgo_North = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,4), S_g_North);

%                     #k_ro = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw)
                    if P_oil(i-1,j) > P_oil(i,j)
                        k_row = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw(i-1,j));
                    else
                        k_row = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw(i,j));
                    end

%                     #k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw)
                    if P_oil(i-1,j) > P_oil(i,j)
                        k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw(i-1,j));
                    else
                        k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw(i,j));
                    end
%                     #println(k_rw)
%                     #k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg)
                    if P_oil(i-1,j) > P_oil(i,j)
                        k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg(i-1,j));
                    else
                        k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg(i,j));
                    end

%                     #k_rog = inearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg)
                    if P_oil(i-1,j) > P_oil(i,j)
                        k_rog = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg(i-1,j));
                    else
                        k_rog = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg(i,j));
                    end
                    k_ro = (k_row + k_rw)*(k_rog + k_rg) - (k_rw + k_rg);
                    B_o = linearInterpolation(oil_prop(:,1), oil_prop(:,3), P_o);
                    mu_o = linearInterpolation(oil_prop(:,1), oil_prop(:,4), P_o);

                    B_w = linearInterpolation(water_prop(:,1), water_prop(:,3), P_o - (P_cow+P_cow_North)/2);
                    mu_w = linearInterpolation(water_prop(:,1), water_prop(:,4), P_o - (P_cow+P_cow_North)/2);

                    B_g = linearInterpolation(gas_prop(:,1), gas_prop(:,3), P_o + (P_cgo+P_cgo_North)/2);
                    mu_g = linearInterpolation(gas_prop(:,1), gas_prop(:,4), P_o + (P_cgo+P_cgo_North)/2);
%                     #first row is oil
                    North_o = calcTran_X(Perm(i-1,j), Perm(i,j), Thickness(i-1,j), Thickness(i,j), DX(i-1,j), DX(i,j), DY(i-1,j), DY(i,j), mu_o, B_o);
                    T(Order(i,j)*3-2, (Order(i-1,j)-1)*3 + 1) = North_o*k_ro;
%                     #second row is gas
                    North_g = calcTran_X(Perm(i-1,j), Perm(i,j), Thickness(i-1,j), Thickness(i,j), DX(i-1,j), DX(i,j), DY(i-1,j), DY(i,j), mu_g, B_g);
                    T(Order(i,j)*3-1, (Order(i-1,j)-1)*3 + 1) = T(Order(i,j)*3-1, (Order(i-1,j)-1)*3 + 1) + North_g*k_rg;
%                     #third row is water
                    North_w = calcTran_X(Perm(i-1,j), Perm(i,j), Thickness(i-1,j), Thickness(i,j), DX(i-1,j), DX(i,j), DY(i-1,j), DY(i,j), mu_w, B_w);
                    T(Order(i,j)*3, (Order(i-1,j)-1)*3 + 1) = North_w*k_rw;


%                     #P_o Sg  and Sw component in matrix of coefficient in 1st row
                    T(Order(i,j)*3-2, Order(i,j)*3-2) = T(Order(i,j)*3-2, Order(i,j)*3-2) - North_o*k_ro;


%                     #P_o Sg  and Sw component in matrix of coefficient in 2nd row
                    R_so_lhs = linearInterpolation(oil_prop(:,1), oil_prop(:,5), P_o);
                    North_g_rso = North_o*k_ro*R_so_lhs;
                    T(Order(i,j)*3-1, (Order(i-1,j)-1)*3 + 1) =  T(Order(i,j)*3-1, (Order(i-1,j)-1)*3 + 1) + North_g_rso;
                    T(Order(i,j)*3-1, Order(i,j)*3-2) = T(Order(i,j)*3-1, Order(i,j)*3-2) - (North_g*k_rg + North_g_rso);


%                     #P_o Sg  and Sw component in matrix of coefficient in 3rd row
                    T(Order(i,j)*3, Order(i,j)*3-2) = T(Order(i,j)*3, Order(i,j)*3-2) - North_w*k_rw;
%                     #gas is missing;

                    Q(Order(i,j)*3-1) = Q(Order(i,j)*3-1) +  (-P_cgo_North + P_cgo)*North_g*k_rg;
                    Q(Order(i,j)*3-0) = Q(Order(i,j)*3-0)  + (P_cow_North - P_cow)*North_w*k_rw;


                end
                if Order(i+1,j) > 0
                    P_o = (P_oil(i+1,j) + P_oil(i,j))/2;
                    S_w_South = Sw(i+1,j);
                    S_g_South = Sg(i+1,j);
                    P_cow_South = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,4), S_w_South);
                    P_cgo_South = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,4), S_g_South);


%                     #k_ro = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw)
                    if P_oil(i+1,j) > P_oil(i,j)
                        k_row = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw(i+1,j));
                    else
                        k_row = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw(i,j));
                    end

%                     #k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw)
                    if P_oil(i+1,j) > P_oil(i,j)
                        k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw(i+1,j));
                    else
                        k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw(i,j));
                    end

%                     #k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg)
                    if P_oil(i+1,j) > P_oil(i,j)
                        k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg(i+1,j));
                    else
                        k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg(i,j));
                    end

%                     #k_rog = inearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg)
                    if P_oil(i+1,j) > P_oil(i,j)
                        k_rog = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg(i+1,j));
                    else
                        k_rog = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg(i,j));
                    end
                    k_ro = (k_row + k_rw)*(k_rog + k_rg) - (k_rw + k_rg);
                    B_o = linearInterpolation(oil_prop(:,1), oil_prop(:,3), P_o);
                    mu_o = linearInterpolation(oil_prop(:,1), oil_prop(:,4), P_o);

                    B_w = linearInterpolation(water_prop(:,1), water_prop(:,3), P_o - (P_cow+P_cow_South)/2);
                    mu_w = linearInterpolation(water_prop(:,1), water_prop(:,4), P_o - (P_cow+P_cow_South)/2);

                    B_g = linearInterpolation(gas_prop(:,1), gas_prop(:,3), P_o + (P_cgo+P_cgo_South)/2);
                    mu_g = linearInterpolation(gas_prop(:,1), gas_prop(:,4), P_o + (P_cgo+P_cgo_South)/2);
                    %first row is oil;
                    South_o = calcTran_X(Perm(i+1,j), Perm(i,j), Thickness(i+1,j), Thickness(i,j), DX(i+1,j), DX(i,j), DY(i+1,j), DY(i,j), mu_o, B_o);
                    T(Order(i,j)*3-2, (Order(i+1,j)-1)*3 + 1) = South_o*k_ro;
                    %second row is gas
                    South_g = calcTran_X(Perm(i+1,j), Perm(i,j), Thickness(i+1,j), Thickness(i,j), DX(i+1,j), DX(i,j), DY(i+1,j), DY(i,j), mu_g, B_g);
                    T(Order(i,j)*3-1, (Order(i+1,j)-1)*3 + 1) = T(Order(i,j)*3-1, (Order(i+1,j)-1)*3 + 1) + South_g*k_rg;
                    %third row is water
                    South_w = calcTran_X(Perm(i+1,j), Perm(i,j), Thickness(i+1,j), Thickness(i,j), DX(i+1,j), DX(i,j), DY(i+1,j), DY(i,j), mu_w, B_w);
                    T(Order(i,j)*3, (Order(i+1,j)-1)*3 + 1) = South_w*k_rw;

                    %P_o Sg  and Sw component in matrix of coefficient in 1st row
                    T(Order(i,j)*3-2, Order(i,j)*3-2) = T(Order(i,j)*3-2, Order(i,j)*3-2) - South_o*k_ro;


                    %P_o Sg  and Sw component in matrix of coefficient in 2nd row
                    R_so_lhs = linearInterpolation(oil_prop(:,1), oil_prop(:,5), P_o);
                    South_g_rso = South_o*k_ro*R_so_lhs;
                    T(Order(i,j)*3-1, (Order(i+1,j)-1)*3 + 1) = T(Order(i,j)*3-1, (Order(i+1,j)-1)*3 + 1) + South_g_rso;
                    T(Order(i,j)*3-1, Order(i,j)*3-2) = T(Order(i,j)*3-1, Order(i,j)*3-2) - (South_g*k_rg + South_g_rso);


                    %P_o Sg  and Sw component in matrix of coefficient in 3rd row
                    T(Order(i,j)*3, Order(i,j)*3-2) = T(Order(i,j)*3, Order(i,j)*3-2)  -South_w*k_rw;
                    %gas is missing
                    Q(Order(i,j)*3-1) = Q(Order(i,j)*3-1) +(-P_cgo_South + P_cgo)*South_g*k_rg;
                    Q(Order(i,j)*3-0) = Q(Order(i,j)*3-0) +(P_cow_South - P_cow)*South_w*k_rw;


                end
                if Order(i,j-1) > 0
                    P_o = (P_oil(i,j-1) + P_oil(i,j))/2;
                    S_w_West = Sw(i,j-1);
                    S_g_West = Sg(i,j-1);
                    P_cow_West = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,4), S_w_West);
                    P_cgo_West = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,4), S_g_West);


                    %k_ro = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw)
                    if P_oil(i,j-1) > P_oil(i,j)
                        k_row = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw(i,j-1));
                    else
                        k_row = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw(i,j));
                    end

                    %k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw)
                    if P_oil(i,j-1) > P_oil(i,j)
                        k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw(i,j-1));
                    else
                        k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw(i,j));
                    end

                    %k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg)
                    if P_oil(i,j-1) > P_oil(i,j)
                        k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg(i,j-1));
                    else
                        k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg(i,j));
                    end

                    %k_rog = inearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg)
                    if P_oil(i,j-1) > P_oil(i,j)
                        k_rog = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg(i,j-1));
                    else
                        k_rog = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg(i,j));
                    end
                    k_ro = (k_row + k_rw)*(k_rog + k_rg) - (k_rw + k_rg);
                    B_o = linearInterpolation(oil_prop(:,1), oil_prop(:,3), P_o);
                    mu_o = linearInterpolation(oil_prop(:,1), oil_prop(:,4), P_o);

                    B_w = linearInterpolation(water_prop(:,1), water_prop(:,3), P_o - (P_cow+P_cow_West)/2);
                    mu_w = linearInterpolation(water_prop(:,1), water_prop(:,4), P_o - (P_cow+P_cow_West)/2);

                    B_g = linearInterpolation(gas_prop(:,1), gas_prop(:,3), P_o + (P_cgo+P_cgo_West)/2);
                    mu_g = linearInterpolation(gas_prop(:,1), gas_prop(:,4), P_o + (P_cgo+P_cgo_West)/2);
                    %first row is oil
                    West_o = calcTran_X(Perm(i,j-1), Perm(i,j), Thickness(i,j-1), Thickness(i,j), DX(i,j-1), DX(i,j), DY(i,j-1), DY(i,j), mu_o, B_o);
                    T(Order(i,j)*3-2, (Order(i,j-1)-1)*3 + 1) = West_o*k_ro;
                    %second row is gas
                    West_g = calcTran_X(Perm(i,j-1), Perm(i,j), Thickness(i,j-1), Thickness(i,j), DX(i,j-1), DX(i,j), DY(i,j-1), DY(i,j), mu_g, B_g);
                    T(Order(i,j)*3-1, (Order(i,j-1)-1)*3 + 1) = T(Order(i,j)*3-1, (Order(i,j-1)-1)*3 + 1) + West_g*k_rg;
                    %third row is water
                    West_w = calcTran_X(Perm(i,j-1), Perm(i,j), Thickness(i,j-1), Thickness(i,j), DX(i,j-1), DX(i,j), DY(i,j-1), DY(i,j), mu_w, B_w);
                    T(Order(i,j)*3, (Order(i,j-1)-1)*3 + 1) = West_w*k_rw;

                    %P_o Sg  and Sw component in matrix of coefficient in 1st row
                    T(Order(i,j)*3-2, Order(i,j)*3-2) = T(Order(i,j)*3-2, Order(i,j)*3-2) - West_o*k_ro;

                    %P_o Sg  and Sw component in matrix of coefficient in 2nd row
                    R_so_lhs = linearInterpolation(oil_prop(:,1), oil_prop(:,5), P_o);
                    West_g_rso = West_o*k_ro*R_so_lhs;
                    T(Order(i,j)*3-1, (Order(i,j-1)-1)*3 + 1) = T(Order(i,j)*3-1, (Order(i,j-1)-1)*3 + 1) + West_g_rso;
                    T(Order(i,j)*3-1, Order(i,j)*3-2) = T(Order(i,j)*3-1, Order(i,j)*3-2) - (West_g*k_rg + West_g_rso);


                    %P_o Sg  and Sw component in matrix of coefficient in 3rd row
                    T(Order(i,j)*3, Order(i,j)*3-2) =T(Order(i,j)*3, Order(i,j)*3-2) - West_w*k_rw;

                    Q(Order(i,j)*3-1) = Q(Order(i,j)*3-1) + (-P_cgo_West + P_cgo)*West_g*k_rg;
                    Q(Order(i,j)*3-0) = Q(Order(i,j)*3-0) + (P_cow_West - P_cow)*West_w*k_rw;

                end
                if Order(i,j+1) > 0
                    P_o = (P_oil(i,j+1) + P_oil(i,j))/2;
                    S_w_East = Sw(i,j+1);
                    S_g_East = Sg(i,j+1);
                    P_cow_East = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,4), S_w_East);
                    P_cgo_East = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,4), S_g_East);


                    %k_ro = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw)
                    if P_oil(i,j+1) > P_oil(i,j)
                        k_row = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw(i,j+1));
                    else
                        k_row = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw(i,j));
                    end

                    %k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw)
                    if P_oil(i,j+1) > P_oil(i,j)
                        k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw(i,j+1));
                    else
                        k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw(i,j));
                    end

                    %k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg)
                    if P_oil(i,j+1) > P_oil(i,j)
                        k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg(i,j+1));
                    else
                        k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg(i,j));
                    end

                    %k_rog = inearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg)
                    if P_oil(i,j+1) > P_oil(i,j)
                        k_rog = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg(i,j+1));
                    else
                        k_rog = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg(i,j));
                    end
                    k_ro = (k_row + k_rw)*(k_rog + k_rg) - (k_rw + k_rg);
                    B_o = linearInterpolation(oil_prop(:,1), oil_prop(:,3), P_o);
                    mu_o = linearInterpolation(oil_prop(:,1), oil_prop(:,4), P_o);

                    B_w = linearInterpolation(water_prop(:,1), water_prop(:,3), P_o - (P_cow+P_cow_East)/2);
                    mu_w = linearInterpolation(water_prop(:,1), water_prop(:,4), P_o - (P_cow+P_cow_East)/2);
                    
                    B_g = linearInterpolation(gas_prop(:,1), gas_prop(:,3), P_o + (P_cgo+P_cgo_East)/2);
                    mu_g = linearInterpolation(gas_prop(:,1), gas_prop(:,4), P_o + (P_cgo+P_cgo_East)/2);
                    %first row is oil
                    East_o = calcTran_X(Perm(i,j+1), Perm(i,j), Thickness(i,j+1), Thickness(i,j), DX(i,j+1), DX(i,j), DY(i,j+1), DY(i,j), mu_o, B_o);
                    T(Order(i,j)*3-2, (Order(i,j+1)-1)*3 + 1) = East_o*k_ro;
                    %second row is gas
                    East_g = calcTran_X(Perm(i,j+1), Perm(i,j), Thickness(i,j+1), Thickness(i,j), DX(i,j+1), DX(i,j), DY(i,j+1), DY(i,j), mu_g, B_g);
                    T(Order(i,j)*3-1, (Order(i,j+1)-1)*3 + 1) = T(Order(i,j)*3-1, (Order(i,j+1)-1)*3 + 1) + East_g*k_rg;
                    %third row is water
                    East_w = calcTran_X(Perm(i,j+1), Perm(i,j), Thickness(i,j+1), Thickness(i,j), DX(i,j+1), DX(i,j), DY(i,j+1), DY(i,j), mu_w, B_w);
                    T(Order(i,j)*3, (Order(i,j+1)-1)*3 + 1) = East_w*k_rw;

                    %P_o Sg  and Sw component in matrix of coefficient in 1st row
                    T(Order(i,j)*3-2, Order(i,j)*3-2) = T(Order(i,j)*3-2, Order(i,j)*3-2) - East_o*k_ro;

                    %P_o Sg  and Sw component in matrix of coefficient in 2nd row
                    R_so_lhs = linearInterpolation(oil_prop(:,1), oil_prop(:,5), P_o);
                    East_g_rso = East_o*k_ro*R_so_lhs;
                    T(Order(i,j)*3-1, (Order(i,j+1)-1)*3 + 1) =T(Order(i,j)*3-1, (Order(i,j+1)-1)*3 + 1) + East_g_rso;
                    T(Order(i,j)*3-1, Order(i,j)*3-2) = T(Order(i,j)*3-1, Order(i,j)*3-2) - (East_g*k_rg + East_g_rso);


                    %P_o Sg  and Sw component in matrix of coefficient in 3rd row
                    T(Order(i,j)*3, Order(i,j)*3-2) = T(Order(i,j)*3, Order(i,j)*3-2) - East_w*k_rw;

                    Q(Order(i,j)*3-1) =  Q(Order(i,j)*3-1) + (-P_cgo_East + P_cgo)*East_g*k_rg;
                    Q(Order(i,j)*3-0) = Q(Order(i,j)*3-0) + (P_cow_East - P_cow)*East_w*k_rw;

                end
                
                

                por_new(i,j) = poro(i,j) + 0*(P_oil (i,j) - P_init(i,j))*poro(i,j);          
                       


                T(Order(i,j)*3-2, Order(i,j)*3-1) = Gamma/B_o_rhs;
                T(Order(i,j)*3-2, Order(i,j)*3-0) = Gamma/B_o_rhs;

                T(Order(i,j)*3-1, Order(i,j)*3-1) = Gamma*(R_so_rhs/B_o_rhs - 1/B_g_rhs);
                T(Order(i,j)*3-1, Order(i,j)*3-0) = Gamma*R_so_rhs/B_o_rhs;

                T(Order(i,j)*3, Order(i,j)*3-0) = -Gamma/B_w_rhs;
               
                B_o_init = linearInterpolation(oil_prop(:,1), oil_prop(:,3), P_init(i,j));
                B_w_init = linearInterpolation(water_prop(:,1), water_prop(:,3), P_init(i,j));
                B_g_init = linearInterpolation(gas_prop(:,1), gas_prop(:,3), P_init(i,j));
                
                R_so_init = linearInterpolation(oil_prop(:,1), oil_prop(:,5), P_init(i,j));

                
                Q(Order(i,j)*3-2) = Q(Order(i,j)*3-2) +  Gamma*(1/B_o_rhs - 1/B_o_init + (Sw_init(i,j) + Sg_init(i,j))/B_o_init);
                Q(Order(i,j)*3-1) = Q(Order(i,j)*3-1) +  Gamma*((R_so_rhs/B_o_rhs - R_so_init/B_o_init) + R_so_init*(Sw_init(i,j) + Sg_init(i,j))/B_o_init) - Gamma*(Sg_init(i,j)/B_g_init);
                Q(Order(i,j)*3-0) = Q(Order(i,j)*3-0) + -Gamma*(Sw_init(i,j)/B_w_init);

            end
        end
    end
    for i = 1:Row
        for j = 1:Col
            if Order(i,j) == 5
                re = 0.198*DX(i,j);
                rw = 0.25;
                P_cow = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,4), Sw(i,j));
                P_cgo = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,4), Sg(i,j));

                B_o_well = linearInterpolation(oil_prop(:,1), oil_prop(:,3), P_oil(i,j));
                mu_o_well = linearInterpolation(oil_prop(:,1), oil_prop(:,4), P_oil(i,j));

                B_w_well = linearInterpolation(water_prop(:,1), water_prop(:,3), P_oil(i,j) - P_cow);
                mu_w_well = linearInterpolation(water_prop(:,1), water_prop(:,4), P_oil(i,j) - P_cow);

                B_g_well = linearInterpolation(gas_prop(:,1), gas_prop(:,3), P_oil(i,j) + P_cgo);
                mu_g_well = linearInterpolation(gas_prop(:,1), gas_prop(:,4), P_oil(i,j) + P_cgo);

                k_row = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,3), Sw(i,j));
                k_rw = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,2), Sw(i,j));
                k_rog = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,3), Sg(i,j));
                k_rg = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,2), Sg(i,j));
                R_so = linearInterpolation(oil_prop(:,1), oil_prop(:,5), P_oil(i,j));
                k_ro_well = (k_row + k_rw)*(k_rog + k_rg) - (k_rw + k_rg);
                P_sf = P_oil(i,j) - 20/(1.127*1e-3*2*pi*Perm(i,j)*50/log(re/rw))/(k_ro_well/B_o_well/mu_o_well);
                Q(Order(i,j)*3 - 2) = Q(Order(i,j)*3 - 2) + 20;
                Q(Order(i,j)*3 - 1) = Q(Order(i,j)*3 - 1) + (1.127*1e-3*2*pi*Perm(i,j)*50/log(re/rw))*(k_rg/(mu_g_well*B_g_well) + R_so*k_ro_well/(mu_o_well*B_o_well))*(P_oil(i,j) - P_sf);
                Q(Order(i,j)*3 - 0) = Q(Order(i,j)*3 - 0) + (1.127*1e-3*2*pi*Perm(i,j)*50/log(re/rw))*(k_rw/(mu_w_well*B_w_well))*(P_oil(i,j) - P_sf);
            end
        end
    end

    out1 = T;
    out2 = Q;
end
