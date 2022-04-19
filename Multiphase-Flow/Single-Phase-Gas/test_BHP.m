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

well_coor = [[5,9]; [9,6]; [4,20]; [8,10]; [9,18]; [6,14]; [11,13]; [7, 22]; [5,17]; [6,5]; [8,14]];
well_spec = [[1, 0]; [1, 0]; [1, 0]; [1, 0]; [1, 0]; [1, 0]; [2,3e6]; [1,0];[1,0];[1,0];[1, 0]];
well_ppt = [[0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0];[0.25,0];[0.25,0] ;[0.25,0];[0.25,0]];
omega = zeros(row,col);
p_wells = zeros(13,length(well_coor));
p_wells(1,:) = ones(length(well_coor),1)*4000;
q_wells = zeros(13,length(well_coor));
q_wells(1,:) = zeros(length(well_coor),1);



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
    for i = 1:length(well_coor)
        ii = well_coor(i,1) + 1;
        jj = well_coor(i,2) + 1;
        p_wells(time+1,i) = Pressure_update(ii,jj);
        q_wells(time+1,i) = omega(ii,jj) * (Pressure_update(ii,jj) - well_spec(i,2));
        Q_total = Q_total + omega(ii,jj) * (Pressure_update(ii,jj) - well_spec(i,2));
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
    set(gcf,'position',[10,10,1000,700]);
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
%     legend('well 1', 'well 2', 'well 3', 'well 4', 'well 5', 'well 6');
    grid on;
end
%legend('well 1', 'well 2', 'well 3', 'well 4', 'well 5', 'well 6');

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



            
            
