clc;
clear;
%% load data
tic;
kx   = load("matlab_cp2_perm.txt");
ky   = load("matlab_cp2_perm.txt");
kx = kx*0.1;
ky = ky*0.1;
dx   = load("matlab_dx.txt");
dy   = load("matlab_dy.txt");
h    = load("matlab_h.txt");
poro = load("matlab_poro.txt");
top = load("matlab_top.txt");

[row, col] = size(kx);
order = zeros(row,col);
cnt= 0;
for i = 1:row
    for j = 1:col
        if kx(i,j) > 0
            cnt = cnt + 1;
            order(i,j) = cnt;
        end
    end
end



%% Initialize pressure and gamma
Pressure = zeros(row,col);
Pressure_init = zeros(row,col);

for i = 1:row
    for j = 1:col
        if order(i,j) > 0
            Pressure(i,j) = 3000;
            Pressure_init(i,j) = 4000;
        end
    end
end
x = zeros(cnt,1);
for i = 1:row
    for j = 1:col
        if order(i,j) > 0
            x(order(i,j)) = Pressure(i,j);
        end
    end
end

well_coor = [[5,9]; [5,11]; [5,13]; [5,15]; [5,17];  [8,9]; [9,11]; [9,13]; [9,15]; [9,17]; [9,19]];
well_coor_ = [[5,9]; [5,11]; [5,13]; [5,15]; [5,17]; [8,9]; [9,11]; [9,13]; [9,15]; [9,17]; [9,19]; [5,21];[5,19]; [5,12];[5,14]; [5,16]; [5,18]; [9,12];[9,14]; [9,16]; [9,18]];
well_coor__ = [[5,9]; [5,11]; [5,13]; [5,15]; [5,17]; [8,9]; [9,11]; [9,13]; [9,15]; [9,17]; [9,19]; [5,21];[5,19]; [7,21]; [5,12];[5,14]; [5,16]; [5,18]; [9,12];[9,14]; [9,16]; [9,18]];
well_coor___ = [[5,9]; [5,11]; [5,13]; [5,15]; [5,17]; [8,9]; [9,11]; [9,13]; [9,15]; [9,17]; [9,19]; [5,21];[5,19]; [7,21];[7,7];  [5,12];[5,14]; [5,16]; [5,18]; [9,12];[9,14]; [9,16]; [9,18]];

well_spec = [[1, 500]; [1, 500]; [1, 500]; [1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500]; [1, 500];[1, 500];];
well_spec_ = [[1, 500]; [1, 500]; [1, 500]; [1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500]; [2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6]];
well_spec__ = [[1, 500]; [1, 500]; [1, 500]; [1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500]; [2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6]];
well_spec___ = [[1, 500]; [1, 500]; [1, 500]; [1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500];[1, 500]; [2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6];[2,10e6]];

well_ppt = [[0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];];
well_ppt_ = [[0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0];[0.25, 0]];
omega = zeros(row,col);
p_wells = zeros(13,length(well_coor__));
p_wells(1,:) = ones(length(well_coor__),1)*4000;
q_wells = zeros(13,length(well_coor__));
q_wells(1,:) = zeros(length(well_coor__),1);

for time = 1:12
    Z_old = gasZ(Pressure_init, order);
    for i = 1:row
        for j = 1:col
            if order(i,j) > 0
                B_old(i,j) = gasFVF(Pressure_init(i,j), Z_old(i,j));
            end
        end
    end
    while true
    [LHS, RHS] = residual(Pressure, Pressure_init, order, cnt, kx, ky, h, dx, dy, poro);
        
        if time <= 8

        for i = 1:length(well_coor)
            ii = well_coor(i,1)+1;
            jj = well_coor(i,2)+1;
            rw = well_ppt(i,1);
            s = well_ppt(i,2);
            re = 0.28*sqrt((ky(ii,jj)/kx(ii,jj))^0.5*dx(ii,jj)^2 + (kx(ii,jj)/ky(ii,jj))^0.5*dy(ii,jj)^2)/...
                ((ky(ii,jj)/kx(ii,jj))^0.25+(kx(ii,jj)/ky(ii,jj))^0.25);
            z = gasZ(Pressure(ii,jj),1);
            rho = gasDensity(Pressure(ii,jj), z);
            miu_well = gasViz(rho);
            B_well = gasFVF(Pressure(ii,jj),z);
            omega(ii,jj) = 2*pi*1.127e-3*sqrt(kx(ii,jj)*ky(ii,jj))*h(ii,jj)/miu_well/B_well/(log(re/rw)+s);
            if well_spec(i,1) == 1
                LHS(order(ii,jj), order(ii,jj)) = LHS(order(ii,jj), order(ii,jj)) - omega(ii,jj);
                RHS(order(ii,jj)) = RHS(order(ii,jj)) -omega(ii,jj)*well_spec(i,2);
            else
                RHS(order(ii,jj)) = RHS(order(ii,jj)) - well_spec(i,2);
            end
        end
        
        elseif time < 11
            for i = 1:length(well_coor_)
            ii = well_coor_(i,1)+1;
            jj = well_coor_(i,2)+1;
            rw = well_ppt_(i,1);
            s = well_ppt_(i,2);
            re = 0.28*sqrt((ky(ii,jj)/kx(ii,jj))^0.5*dx(ii,jj)^2 + (kx(ii,jj)/ky(ii,jj))^0.5*dy(ii,jj)^2)/...
                ((ky(ii,jj)/kx(ii,jj))^0.25+(kx(ii,jj)/ky(ii,jj))^0.25);
            z = gasZ(Pressure(ii,jj),1);
            rho = gasDensity(Pressure(ii,jj), z);
            miu_well = gasViz(rho);
            B_well = gasFVF(Pressure(ii,jj),z);
            omega(ii,jj) = 2*pi*1.127e-3*sqrt(kx(ii,jj)*ky(ii,jj))*h(ii,jj)/miu_well/B_well/(log(re/rw)+s);
            if well_spec_(i,1) == 1
                LHS(order(ii,jj), order(ii,jj)) = LHS(order(ii,jj), order(ii,jj)) - omega(ii,jj);
                RHS(order(ii,jj)) = RHS(order(ii,jj)) -omega(ii,jj)*well_spec_(i,2);
            else
                RHS(order(ii,jj)) = RHS(order(ii,jj)) - well_spec_(i,2);
            end
            end
        elseif time < 12
            for i = 1:length(well_coor__)
            ii = well_coor__(i,1)+1;
            jj = well_coor__(i,2)+1;
            rw = well_ppt_(i,1);
            s = well_ppt_(i,2);
            re = 0.28*sqrt((ky(ii,jj)/kx(ii,jj))^0.5*dx(ii,jj)^2 + (kx(ii,jj)/ky(ii,jj))^0.5*dy(ii,jj)^2)/...
                ((ky(ii,jj)/kx(ii,jj))^0.25+(kx(ii,jj)/ky(ii,jj))^0.25);
            z = gasZ(Pressure(ii,jj),1);
            rho = gasDensity(Pressure(ii,jj), z);
            miu_well = gasViz(rho);
            B_well = gasFVF(Pressure(ii,jj),z);
            omega(ii,jj) = 2*pi*1.127e-3*sqrt(kx(ii,jj)*ky(ii,jj))*h(ii,jj)/miu_well/B_well/(log(re/rw)+s);
            if well_spec__(i,1) == 1
                LHS(order(ii,jj), order(ii,jj)) = LHS(order(ii,jj), order(ii,jj)) - omega(ii,jj);
                RHS(order(ii,jj)) = RHS(order(ii,jj)) -omega(ii,jj)*well_spec__(i,2);
            else
                RHS(order(ii,jj)) = RHS(order(ii,jj)) - well_spec__(i,2);
            end
            end
        else
            for i = 1:length(well_coor___)
            ii = well_coor___(i,1)+1;
            jj = well_coor___(i,2)+1;
            rw = well_ppt_(i,1);
            s = well_ppt_(i,2);
            re = 0.28*sqrt((ky(ii,jj)/kx(ii,jj))^0.5*dx(ii,jj)^2 + (kx(ii,jj)/ky(ii,jj))^0.5*dy(ii,jj)^2)/...
                ((ky(ii,jj)/kx(ii,jj))^0.25+(kx(ii,jj)/ky(ii,jj))^0.25);
            z = gasZ(Pressure(ii,jj),1);
            rho = gasDensity(Pressure(ii,jj), z);
            miu_well = gasViz(rho);
            B_well = gasFVF(Pressure(ii,jj),z);
            omega(ii,jj) = 2*pi*1.127e-3*sqrt(kx(ii,jj)*ky(ii,jj))*h(ii,jj)/miu_well/B_well/(log(re/rw)+s);
            if well_spec___(i,1) == 1
                LHS(order(ii,jj), order(ii,jj)) = LHS(order(ii,jj), order(ii,jj)) - omega(ii,jj);
                RHS(order(ii,jj)) = RHS(order(ii,jj)) -omega(ii,jj)*well_spec___(i,2);
            else
                RHS(order(ii,jj)) = RHS(order(ii,jj)) - well_spec___(i,2);
            end
            end
            
            
        end
        
            
        Residual = LHS*x - RHS;


        J = zeros(cnt,cnt);
        for i = 1:row
            for j = 1:col
                if order(i,j) > 0 
                    P_new = Pressure;
                    P_new(i,j) = P_new(i,j) + 1e-6;
                    x_new = x;
                    x_new(order(i,j)) = x_new(order(i,j))+1e-6;
                    [LHS_new, RHS_new] = residual(P_new, Pressure_init, order, cnt, kx, ky, h, dx, dy, poro);
                if time <= 8
                for k = 1:length(well_coor)
                    ii = well_coor(k,1)+1;
                    jj = well_coor(k,2)+1;
                    if well_spec(k,1) == 1
                        LHS_new(order(ii,jj), order(ii,jj)) = LHS_new(order(ii,jj), order(ii,jj)) - omega(ii,jj);
                        RHS_new(order(ii,jj)) = RHS_new(order(ii,jj)) -omega(ii,jj)*well_spec(k,2);
                    else
                        RHS_new(order(ii,jj)) = RHS_new(order(ii,jj)) - well_spec(k,2);
                    end
                end
                elseif time < 11
                    for k = 1:length(well_coor_)
                    ii = well_coor_(k,1)+1;
                    jj = well_coor_(k,2)+1;
                    if well_spec_(k,1) == 1
                        LHS_new(order(ii,jj), order(ii,jj)) = LHS_new(order(ii,jj), order(ii,jj)) - omega(ii,jj);
                        RHS_new(order(ii,jj)) = RHS_new(order(ii,jj)) -omega(ii,jj)*well_spec_(k,2);
                    else
                        RHS_new(order(ii,jj)) = RHS_new(order(ii,jj)) - well_spec_(k,2);
                    end
                    end
                elseif time < 12
                    for k = 1:length(well_coor__)
                    ii = well_coor__(k,1)+1;
                    jj = well_coor__(k,2)+1;
                    if well_spec__(k,1) == 1
                        LHS_new(order(ii,jj), order(ii,jj)) = LHS_new(order(ii,jj), order(ii,jj)) - omega(ii,jj);
                        RHS_new(order(ii,jj)) = RHS_new(order(ii,jj)) -omega(ii,jj)*well_spec__(k,2);
                    else
                        RHS_new(order(ii,jj)) = RHS_new(order(ii,jj)) - well_spec__(k,2);
                    end
                    end
                else
                    for k = 1:length(well_coor___)
                    ii = well_coor___(k,1)+1;
                    jj = well_coor___(k,2)+1;
                    if well_spec___(k,1) == 1
                        LHS_new(order(ii,jj), order(ii,jj)) = LHS_new(order(ii,jj), order(ii,jj)) - omega(ii,jj);
                        RHS_new(order(ii,jj)) = RHS_new(order(ii,jj)) -omega(ii,jj)*well_spec___(k,2);
                    else
                        RHS_new(order(ii,jj)) = RHS_new(order(ii,jj)) - well_spec___(k,2);
                    end
                    end
                    
                end
                    
                
                Residual_new = LHS_new*x_new - RHS_new;
                J(:,order(i,j)) = (Residual_new - Residual)/1e-6;
                end
            end
        end
    
        result = J\Residual;
        if max(abs(result)) < 0.01
            break
        end

        disp(result(1));
        Pressure_update = zeros(row,col);
        for i = 1:row
            for j = 1:col
                if order(i,j) > 0
                    Pressure_update(i,j) = Pressure(i,j) - result(order(i,j));
                end
            end
        end

        Pressure = Pressure_update;

    
        for i = 1:row
            for j = 1:col
                if order(i,j) > 0
                    x(order(i,j)) = Pressure(i,j);
                end
            end
        end
    end
    
    Z_new = gasZ(Pressure_update, order);
    for i = 1:row
        for j = 1:col
            if order(i,j) > 0
                B_new(i,j) = gasFVF(Pressure_update(i,j), Z_new(i,j));
            end
        end
    end
    IMBC = 0;
    Q_total = 0;
    if time <= 8
    for i = 1:length(well_coor)
        ii = well_coor(i,1) + 1;
        jj = well_coor(i,2) + 1;
        if well_spec(i,1) == 1
            p_wells(time+1,i) = Pressure_update(ii,jj);
            q_wells(time+1,i) = omega(ii,jj) * (Pressure_update(ii,jj) - well_spec(i,2));
            Q_total = Q_total + omega(ii,jj) * (Pressure_update(ii,jj) - well_spec(i,2));
        else
            p_wells(time+1,i) = Pressure_update(ii,jj);
            q_wells(time+1,i) =  - well_spec(i,2);
            Q_total = Q_total - well_spec(i,2);
        end      
    end
    elseif time < 11
        for i = 1:length(well_coor_)
        ii = well_coor_(i,1) + 1;
        jj = well_coor_(i,2) + 1;
        if well_spec_(i,1) == 1
            p_wells(time+1,i) = Pressure_update(ii,jj);
            q_wells(time+1,i) = omega(ii,jj) * (Pressure_update(ii,jj) - well_spec_(i,2));
            Q_total = Q_total + omega(ii,jj) * (Pressure_update(ii,jj) - well_spec_(i,2));
        else
            p_wells(time+1,i) = Pressure_update(ii,jj);
            q_wells(time+1,i) =  - well_spec_(i,2);
            Q_total = Q_total - well_spec_(i,2);
        end      
        end
    elseif time < 12
        for i = 1:length(well_coor__)
        ii = well_coor__(i,1) + 1;
        jj = well_coor__(i,2) + 1;
        if well_spec__(i,1) == 1
            p_wells(time+1,i) = Pressure_update(ii,jj);
            q_wells(time+1,i) = omega(ii,jj) * (Pressure_update(ii,jj) - well_spec__(i,2));
            Q_total = Q_total + omega(ii,jj) * (Pressure_update(ii,jj) - well_spec__(i,2));
        else
            p_wells(time+1,i) = Pressure_update(ii,jj);
            q_wells(time+1,i) =  - well_spec__(i,2);
            Q_total = Q_total - well_spec__(i,2);
        end      
        end
    else
        for i = 1:length(well_coor___)
        ii = well_coor___(i,1) + 1;
        jj = well_coor___(i,2) + 1;
        if well_spec___(i,1) == 1
            p_wells(time+1,i) = Pressure_update(ii,jj);
            q_wells(time+1,i) = omega(ii,jj) * (Pressure_update(ii,jj) - well_spec___(i,2));
            Q_total = Q_total + omega(ii,jj) * (Pressure_update(ii,jj) - well_spec___(i,2));
        else
            p_wells(time+1,i) = Pressure_update(ii,jj);
            q_wells(time+1,i) =  - well_spec___(i,2);
            Q_total = Q_total - well_spec___(i,2);
        end      
        end
        
    end
    
 

    for i =1:row
        for j = 1:col
            if order(i,j) > 0
                IMBC = IMBC + dx(i,j)*dy(i,j)*h(i,j)*poro(i,j)/5.615*(1/B_new(i,j) - 1/B_old(i,j));
            end
        end
    end
    fprintf("IMBC check %0.8f \n",abs(IMBC/Q_total/30));
    figure(1);
    heatmap(Pressure_update);
    %pcolor(Pressure_update);
    colormap jet;
    %s.LineWidth = 0.5;
    %axis square;
    %imagesc(Pressure_update);
    colorbar;
    str = sprintf('%d',time);
    set(gcf,'position',[10,10,1000,700]);
    title("Block pressure distribution at" + " " + str + " " +"month");
    pause(0.2);
    

    Pressure_init = Pressure;

end

%% Plotting pressure and production profile for all wells.
time = linspace(0,12,13);
%line_color = ['b' 'g' 'r' 'c' 'm' 'y'];


for i = 1:length(well_coor)
    figure(2);
    title("Pressure history of all wells");
    ylabel("Pressure (psi)");
    xlabel("Times (moths)");
    xlim([0, max(time)+1]);
    ylim([1000, 5000]);
    hold on;
    plot(time, p_wells(:,i), 'LineWidth',3);
    grid on;
end
%legend('well 1', 'well 2', 'well 3', 'well 4', 'well 5', 'well 6');
Q_wells = q_wells.*(q_wells > 0);
Q_total_monthly = sum(Q_wells,2);
for i = 1:length(well_coor)
    figure(3);
    title("Production history of all wells");
    ylabel("Production rate (SCF/D)");
    xlabel("Times (months)");
    xlim([0, max(time)+1]);
    ylim([0, 30e6]);
    hold on;
    plot(time, q_wells(:,i), 'LineWidth',3);
    grid on;
end
%legend('well 1', 'well 2', 'well 3', 'well 4', 'well 5', 'well 6');
figure(4);
plot(time, Q_total_monthly*30, "o-", 'LineWidth', 3,'markersize',10,'markerfacecolor','g');
title("Total production of each month in a year");
ylabel("Production (SCF)");
xlabel("Times (months)");

demand = 1e3*[1   140   130   115    90    60    40    28    25    45    65    85   120]';
Number_households = Q_total_monthly*30./demand;
figure(5);
plot(time, Number_households, "o-", 'LineWidth', 3,'markersize',10,'markerfacecolor','g');
title("Total number of households meet the gas demand");
ylabel("Number of households");
xlabel("Times (months)");

            
            
