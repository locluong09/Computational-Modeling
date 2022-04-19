clc;
clear;
format long;

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

h = boundaries(Thick);
dx = boundaries(DX);
dy = boundaries(DY);
poro = boundaries(Phi);
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
iter = 0;

% for i = 1:1
while true
    iter = iter + 1;
    [T, Q]  = construct_Residual_matrix(P_oil, Sg, Sw, Order, Pre_init, Sg_init, Sw_init);
    Residual = T*x - Q;
    J = Jacobian(P_oil, Sg, Sw, Order, Pre_init, Sg_init, Sw_init);
    Delta = J\Residual;

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
    if max(abs(Delta)) < 0.0001
         break
    end
    %if maximum(broadcast(abs, x_new-x)) < 0.01
        %break
    %end

    P_oil = P_update;
    Sg = Sg_update;
    Sw = Sw_update;
    x = x_new;

end


for i = 1:Row
    for j = 1:Col
        if Order(i,j) > 0
            B_oil(i,j) = linearInterpolation(oil_prop(:,1), oil_prop(:,3), P_oil(i,j));
            P_cow_ = linearInterpolation(oil_water_rel(:,1), oil_water_rel(:,4), Sw(i,j));
            P_cgo_ = linearInterpolation(gas_oil_rel(:,1), gas_oil_rel(:,4), Sg(i,j));
            B_gas(i,j) = linearInterpolation(gas_prop(:,1), gas_prop(:,3), P_oil(i,j) + P_cgo_);
            B_water(i,j) = linearInterpolation(water_prop(:,1), water_prop(:,3), P_oil(i,j) - P_cow_);
            B_oil_init(i,j) = linearInterpolation(oil_prop(:,1), oil_prop(:,3), 4800);
            B_gas_init(i,j) = linearInterpolation(gas_prop(:,1), gas_prop(:,3), 4800);
            B_water_init(i,j) = linearInterpolation(water_prop(:,1), water_prop(:,3), 4800);
            R_so_mb(i,j) = linearInterpolation(oil_prop(:,1), oil_prop(:,5), P_oil(i,j));
            R_so_old(i,j) = linearInterpolation(oil_prop(:,1), oil_prop(:,5), 4800);
        end
    end
end

IMBC_oil = 0;
IMBC_gas = 0;
IMBC_water = 0;

for i =1:Row
    for j = 1:Col
        if Order(i,j) > 0
            IMBC_oil = IMBC_oil + dx(i,j)*dy(i,j)*h(i,j)*0.2/5.615*((1-Sg(i,j)-Sw(i,j))/B_oil(i,j) - (1-Sw_init(i,j)-Sg_init(i,j))/B_oil_init(i,j));
            IMBC_gas = IMBC_gas + dx(i,j)*dy(i,j)*h(i,j)*0.2/5.615*(Sg(i,j)/B_gas(i,j) - Sg_init(i,j)/B_gas_init(i,j)) + dx(i,j)*dy(i,j)*h(i,j)*0.2/5.615*(R_so_mb(i,j)*(1-Sg(i,j)-Sw(i,j))/B_oil(i,j) - R_so_old(i,j)*(1-Sw_init(i,j)-Sg_init(i,j))/B_oil_init(i,j));
            IMBC_water = IMBC_water + dx(i,j)*dy(i,j)*h(i,j)*0.2/5.615*(Sw(i,j)/B_water(i,j) - Sw_init(i,j)/B_water_init(i,j));
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
                Q_oil = 20;
                Q_gas = (1.127*1e-3*2*pi*Perm(i,j)*50/log(re/rw))*(k_rg/(mu_g_well*B_g_well) + R_so*k_ro_well/(mu_o_well*B_o_well))*(P_oil(i,j) - P_sf);
                Q_water =  (1.127*1e-3*2*pi*Perm(i,j)*50/log(re/rw))*(k_rw/(mu_w_well*B_w_well))*(P_oil(i,j) - P_sf);
        end
    end
end

fprintf("IMBC oil check %0.8f \n",abs(IMBC_oil/Q_oil/5));

fprintf("IMBC gas check %0.8f \n",abs(IMBC_gas/Q_gas/5));

fprintf("IMBC water check %0.8f \n",abs(IMBC_water/Q_water/5));


