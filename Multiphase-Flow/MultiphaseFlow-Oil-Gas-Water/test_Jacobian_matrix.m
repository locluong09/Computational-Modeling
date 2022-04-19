clc;
clear;
format long;



Perm = ones(3,3)*100;
Perm = boundaries(Perm);
[Row, Col] = size(Perm);
cnt = 0;
Order = zeros(Row, Col);
for i = 1:Row
    for j = 1:Col
        if Perm(i,j) > 0
            cnt = cnt+ 1;
            Order(i,j) = cnt;
        end
    end
end

Pre_init = ones(3,3)*4800;
Sw_init = ones(3,3)*0.3;
Sg_init = ones(3,3)*0.1;

Pre_init = boundaries(Pre_init);
Sw_init = boundaries(Sw_init);
Sg_init = boundaries(Sg_init);
% initial pressure of oil phase
P_oil = boundaries(ones(3,3)*4800);
Sw = boundaries(ones(3,3)*0.3);
Sg = boundaries(ones(3,3)*0.1);

x = zeros(27,1);

for i = 1:Row
    for j = 1:Col
        if Order(i,j) > 0
            x(Order(i,j)*3-2) = P_oil(i,j);
            x(Order(i,j)*3-1) = Sg(i,j);
            x(Order(i,j)*3) = Sw(i,j);
        end
    end
end


for i = 1:3
    [T, Q]  = construct_Residual_matrix(P_oil, Sg, Sw, Order, Pre_init, Sg_init, Sw_init);
    Residual = T*x - Q;
    Jacobian_matrix = Jacobian(P_oil, Sg, Sw, Order, Pre_init, Sg_init, Sw_init);
    Delta = Jacobian_matrix\Residual;

    x_new = x;
    
    x_new = x_new - Delta;
    
   
    P_update = zeros(Row,Col);
    Sg_update = zeros(Row,Col);
    Sw_update = zeros(Row,Col);
    
    for i = 1:Row
        for j = 1:Col
            if Order(i,j) > 0
                P_update(i,j) = x_new(Order(i,j)*3-2);
                Sg_update(i,j) = x_new(Order(i,j)*3-1);
                Sw_update(i,j) = x_new(Order(i,j)*3);
            end
        end
    end

    P_oil = P_update;
    Sg = Sg_update;
    Sw = Sw_update;
    x = x_new;

end



