clc;
clear;
%% load data
tic;
kx   = load("cpp_kx.txt");
ky   = load("cpp_ky.txt");
dx   = load("cpp_dx.txt");
dy   = load("cpp_dy.txt");
h    = load("cpp_h.txt");
poro = load("cpp_poro.txt");
top = load("cpp_top.txt");

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

%% Initialize matrix N,S,W,E
N = zeros(row, col);
S = zeros(row, col);
W = zeros(row, col);
E = zeros(row, col);
N_ = zeros(row, col);
S_ = zeros(row, col);
W_ = zeros(row, col);
E_ = zeros(row, col);
C = zeros(row,col);
Q = zeros(row,col);

%% Initialize pressure and gamma
Pressure = zeros(row,col);
for i = 1:row
    for j = 1:col
        if order(i,j) > 0
            Pressure(i,j) = 5000;
        end
    end
end
gamma = zeros(row, col);
B_0 = arrayfun(@(Pressure) B_avg(Pressure), Pressure);
Por_0 = poro;
%% Save pressure block and production rate

p_wells = zeros(36,6);
q_wells = zeros(36,6);
mbal1 = zeros(36,1);
mbal2 = zeros(36,1);
%% Construct coefficient matrix and solving the pressure distribution
for t = 1:36
    for i = 1:row
        for j = 1:col
            if order(i,j) > 0
                if order(i-1,j) > 0
                    P_avg = (Pressure(i,j) + Pressure(i-1,j))/2;
                    mu = mu_avg(P_avg);
                    B = B_avg(P_avg);
                    rho = rho_avg(P_avg);
                    N(i,j) = tran_X(kx(i,j), kx(i-1,j), h(i,j), h(i-1,j), mu, B, dx(i,j), dx(i-1,j), dy(i,j), dy(i-1,j));
                    N_(i,j) = N(i,j)*rho/144;
                end
                
                if order(i+1,j) > 0
                    P_avg = (Pressure(i,j) + Pressure(i+1,j))/2;
                    mu = mu_avg(P_avg);
                    B = B_avg(P_avg);
                    rho = rho_avg(P_avg);
                    S(i,j) = tran_X(kx(i,j), kx(i+1,j), h(i,j), h(i+1,j), mu, B, dx(i,j), dx(i+1,j), dy(i,j), dy(i+1,j));
                    S_(i,j) = S(i,j)*rho/144;
                end
                
                if order(i,j-1) > 0
                    P_avg = (Pressure(i,j) + Pressure(i,j-1))/2;
                    mu = mu_avg(P_avg);
                    B = B_avg(P_avg);
                    rho = rho_avg(P_avg);
                    W(i,j) = tran_Y(ky(i,j), ky(i,j-1), h(i,j), h(i,j-1), mu, B, dx(i,j), dx(i,j-1), dy(i,j), dy(i,j-1));
                    W_(i,j) = W(i,j)*rho/144;
                end
                
                if order(i,j+1) > 0
                    P_avg = (Pressure(i,j) + Pressure(i,j+1))/2;
                    mu = mu_avg(P_avg);
                    B = B_avg(P_avg);
                    rho = rho_avg(P_avg);
                    E(i,j) = tran_Y(ky(i,j), ky(i,j+1), h(i,j), h(i,j+1), mu, B, dx(i,j), dx(i,j+1), dy(i,j), dy(i,j+1));
                    E_(i,j) = E(i,j)*rho/144;
                end
                gamma(i,j) = dx(i,j)*dy(i,j)*h(i,j)*poro(i,j)*3.1e-5/(5.615*30);
                C(i,j) = -(N(i,j) + S(i,j) + W(i,j) + E(i,j) +gamma(i,j));
                Q(i,j) = N_(i,j)*(top(i-1,j) - top(i,j)) + S_(i,j)*(top(i+1,j) - top(i,j))+W_(i,j)*(top(i,j-1) - top(i,j)) + E_(i,j)*(top(i,j+1) - top(i,j)) - gamma(i,j)*Pressure(i,j);
            end
        end
    end
%% Determine FVF of oil
    B_init = arrayfun(@(Pressure) B_avg(Pressure), Pressure);
    
                
%% Well specification

    well_coor = [[4,6]; [9,6]; [4,20]; [8,10]; [12,13]; [6,14]];
    well_spec = [[1, 1000]; [1, 1000]; [1, 1000]; [1, 1000]; [1, 1000]; [1, 1000]];
    well_ppt = [[0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]];
    omega = zeros(row,col);

    if t < 24
        for i = 1:length(well_coor)
            ii = well_coor(i,1)+1;
            jj = well_coor(i,2)+1;
            rw = well_ppt(i,1);
            s = well_ppt(i,2);
            re = 0.28*sqrt((ky(ii,jj)/kx(ii,jj))^0.5*dx(ii,jj)^2 + (kx(ii,jj)/ky(ii,jj))^0.5*dy(ii,jj)^2)/...
                ((ky(ii,jj)/kx(ii,jj))^0.25+(kx(ii,jj)/ky(ii,jj))^0.25);
            miu_well = mu_avg(Pressure(ii,jj));
            B_well = B_avg(Pressure(ii,jj));
            omega(ii,jj) = 2*pi*1.127e-3*sqrt(kx(ii,jj)*ky(ii,jj))*h(ii,jj)/miu_well/B_well/(log(re/rw)+s);
            if well_spec(i,1) == 1
                C(ii,jj) = C(ii,jj) - omega(ii,jj);
                Q(ii,jj) = Q(ii,jj) -omega(ii,jj)*well_spec(i,2);
            else
                Q(ii,jj) = Q(ii,jj) - well_spec(i,2);
            end
            
        end
    else
        for i = 3:length(well_coor)
            ii = well_coor(i,1)+1;
            jj = well_coor(i,2)+1;
            rw = well_ppt(i,1);
            s = well_ppt(i,2);
            re = 0.28*sqrt((ky(ii,jj)/kx(ii,jj))^0.5*dx(ii,jj)^2 + (kx(ii,jj)/ky(ii,jj))^0.5*dy(ii,jj)^2)/...
                ((ky(ii,jj)/kx(ii,jj))^0.25+(kx(ii,jj)/ky(ii,jj))^0.25);
            miu_well = mu_avg(Pressure(ii,jj));
            B_well = B_avg(Pressure(ii,jj));
            omega(ii,jj) = 2*pi*1.127e-3*sqrt(kx(ii,jj)*ky(ii,jj))*h(ii,jj)/miu_well/B_well/(log(re/rw)+s);
            if well_spec(i,1) == 1
                C(ii,jj) = C(ii,jj) - omega(ii,jj);
                Q(ii,jj) = Q(ii,jj) -omega(ii,jj)*well_spec(i,2);
            else
                Q(ii,jj) = Q(ii,jj) - well_spec(i,2);
            end
        end
    end
    
        


%% Construct LHS and RHS
    LHS = zeros(cnt, cnt);
    RHS = zeros(cnt,1);
    for i =1:row
        for j = 1:col
            if order(i,j)>0
                LHS(order(i,j), order(i,j)) = C(i,j);
                if order(i-1,j) > 0
                    LHS(order(i,j), order(i-1,j)) = N(i,j);
                end
                
                if order(i+1,j) > 0
                    LHS(order(i,j), order(i+1,j)) = S(i,j);
                end
                
                if order(i,j-1) > 0
                    LHS(order(i,j), order(i,j-1)) = W(i,j);
                end
                
                if order(i,j+1) > 0
                    LHS(order(i,j), order(i,j+1)) = E(i,j);
                end
                RHS(order(i,j)) = Q(i,j);
            end
        end
    end

%% Update pressure
    x = LHS\RHS;
    Pressure_update = zeros(row,col);
    por_new = zeros(row,col);

    for i = 1:row
        for j = 1:col
            if order(i,j) > 0
                Pressure_update(i,j) = x(order(i,j));
                por_new(i,j) = poro(i,j) + 1e-6*(Pressure_update (i,j) - Pressure(i,j))*poro(i,j);          
            end
        end
    end
%% Plot pressure update
%{
    figure(1);
    heatmap(Pressure_update);
    colormap jet;
    %s.LineWidth = 0.5;
    %axis square;
    %imagesc(Pressure_update);
    colorbar;
    str = sprintf('%d',30*t);
    title("Block pressure distribution at" + " " + str + " " +"days");
    pause(0.5);   
%}    
    Q_total = 0;
    for i = 1:length(well_coor)
        ii = well_coor(i,1) + 1;
        jj = well_coor(i,2) + 1;
        p_wells(t,i) = Pressure_update(ii,jj);
        q_wells(t,i) = omega(ii,jj) * (Pressure_update(ii,jj) - 1000);
        Q_total = Q_total + omega(ii,jj) * (Pressure_update(ii,jj) - 1000);
    end
       
    

%% Material balance check            
    B_new = arrayfun(@(Pressure_update) B_avg(Pressure_update), Pressure_update);
    IMBC = 0;
    for i =1:row
        for j = 1:col
            if order(i,j) > 0
                IMBC = IMBC + dx(i,j)*dy(i,j)*h(i,j)/5.615*(por_new(i,j)/B_new(i,j) - poro(i,j)/B_init(i,j));
            end
        end
    end
    fprintf("IMBC check %0.8f \n",abs(IMBC/Q_total/30));
    mbal1(t) = abs(IMBC/Q_total/30);
    
    CMBC = 0;
    for i =1:row
        for j = 1:col
            if order(i,j) > 0
                CMBC = CMBC + dx(i,j)*dy(i,j)*h(i,j)/5.615*(por_new(i,j)/B_new(i,j) - Por_0(i,j)/B_0(i,j));
            end
        end
    end
    total_Q = sum(sum(q_wells));
    fprintf("CMBC check %0.8f \n",abs(CMBC/total_Q/30));
    mbal2(t) = abs(CMBC/total_Q/30);
    Pressure = Pressure_update;
    poro = por_new;
    
end

%% Plotting pressure and production profile for all wells.
%{
time = linspace(30,360*3,36);
line_color = ['b' 'g' 'r' 'c' 'm' 'y'];


for i = 1:6
    figure(2);
    title("Pressure history of all wells");
    ylabel("Pressure (psi)");
    xlabel("Times (days)");
    xlim([0, max(time) + 30]);
    ylim([1000, 5000]);
    hold on;
    plot(time, p_wells(:,i), line_color(i), 'LineWidth',3);
%     legend('well 1', 'well 2', 'well 3', 'well 4', 'well 5', 'well 6');
    grid on;
end
legend('well 1', 'well 2', 'well 3', 'well 4', 'well 5', 'well 6');

for i = 1:6
    figure(3);
    title("Production history of all wells");
    ylabel("Production rate (STB/D)");
    xlabel("Times (days)");
    xlim([0, max(time) + 30]);
    hold on;
    plot(time, q_wells(:,i), line_color(i), 'LineWidth',3);
    grid on;
end
legend('well 1', 'well 2', 'well 3', 'well 4', 'well 5', 'well 6');
            
figure(4)
plot(time, mbal1, 'r');
title("Incremental material balance check");
ylabel("IMBC");
xlabel("Time");
grid on;

figure(5)
plot(time, mbal2, 'b');
title("Cumulative material balance check");
ylabel("CMBC");
xlabel("Time");
grid on;
%}

runtime1 = toc;
fprintf("Time taken: %0.6f", runtime1)





