function [out1, out2] = residual(Pressure,Pressure_old, order, cnt, kx, ky, h, dx, dy, poro)
    [row,col] = size(order);
    N = zeros(row, col);
    S = zeros(row, col);
    W = zeros(row, col);
    E = zeros(row, col);
   

    C = zeros(row,col);
    Q = zeros(row,col);
    Z = gasZ(Pressure, order);
    
    
    for i = 1:row
        for j = 1:col
            if order(i,j) > 0
                if order(i-1,j) > 0
                    P_avg = (Pressure(i,j) + Pressure(i-1,j))/2;
                    z = gasZ(P_avg,1);
					rho = gasDensity(P_avg, z);
					mu = gasViz(rho);
					B = gasFVF(P_avg, z);
                    N(i,j) = tran_X(kx(i,j), kx(i-1,j), h(i,j), h(i-1,j), mu, B, dx(i,j), dx(i-1,j), dy(i,j), dy(i-1,j));
                    
                end
                
                if order(i+1,j) > 0
                    P_avg = (Pressure(i,j) + Pressure(i+1,j))/2;
                    z = gasZ(P_avg,1);
					rho = gasDensity(P_avg, z);
					mu = gasViz(rho);
					B = gasFVF(P_avg, z);
                    S(i,j) = tran_X(kx(i,j), kx(i+1,j), h(i,j), h(i+1,j), mu, B, dx(i,j), dx(i+1,j), dy(i,j), dy(i+1,j));
                    
                end
                
                if order(i,j-1) > 0
                    P_avg = (Pressure(i,j) + Pressure(i,j-1))/2;
                    z = gasZ(P_avg,1);
					rho = gasDensity(P_avg, z);
					mu = gasViz(rho);
					B = gasFVF(P_avg, z);
                    W(i,j) = tran_Y(ky(i,j), ky(i,j-1), h(i,j), h(i,j-1), mu, B, dx(i,j), dx(i,j-1), dy(i,j), dy(i,j-1));
                    
                end
                
                if order(i,j+1) > 0
                    P_avg = (Pressure(i,j) + Pressure(i,j+1))/2;
                    z = gasZ(P_avg,1);
					rho = gasDensity(P_avg, z);
					mu = gasViz(rho);
					B = gasFVF(P_avg, z);
                    E(i,j) = tran_Y(ky(i,j), ky(i,j+1), h(i,j), h(i,j+1), mu, B, dx(i,j), dx(i,j+1), dy(i,j), dy(i,j+1));
                    
                end
                gamma(i,j) = dx(i,j)*dy(i,j)*h(i,j)*poro(i,j)*520/(14.7*(150+460)*30);
                C(i,j) = -(N(i,j) + S(i,j) + W(i,j) + E(i,j) +gamma(i,j)/Z(i,j));
                Q(i,j) = gamma(i,j)*(-Pressure_old(i,j)/gasZ(Pressure_old(i,j),1));
            end
        end
    end
    
    %Q(9,11) = Q(9,11) + 30e6;
    
%     well_coor = [[4,6]; [9,6]; [4,20]; [8,10]; [12,13]; [6,14]];
%     well_spec = [[2, -3e6]; [2, -3e6]; [2, -3e6]; [2, -3e6]; [2, -3e6]; [2, -3e6]];
%     well_ppt = [[0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]; [0.25, 0]];
%     omega = zeros(row,col);
%     
%     for i = 1:length(well_coor)
%             ii = well_coor(i,1)+1;
%             jj = well_coor(i,2)+1;
%             rw = well_ppt(i,1);
%             s = well_ppt(i,2);
%             re = 0.28*sqrt((ky(ii,jj)/kx(ii,jj))^0.5*dx(ii,jj)^2 + (kx(ii,jj)/ky(ii,jj))^0.5*dy(ii,jj)^2)/...
%                 ((ky(ii,jj)/kx(ii,jj))^0.25+(kx(ii,jj)/ky(ii,jj))^0.25);
%             miu_well = mu_avg(Pressure(ii,jj));
%             B_well = B_avg(Pressure(ii,jj));
%             omega(ii,jj) = 2*pi*1.127e-3*sqrt(kx(ii,jj)*ky(ii,jj))*h(ii,jj)/miu_well/B_well/(log(re/rw)+s);
%             if well_spec(i,1) == 1
%                 C(ii,jj) = C(ii,jj) - omega(ii,jj);
%                 Q(ii,jj) = Q(ii,jj) -omega(ii,jj)*well_spec(i,2);
%             else
%                 Q(ii,jj) = Q(ii,jj) - well_spec(i,2);
%             end
%     end
    
    
    
    
    
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
    out1 = LHS;
    out2 = RHS;
end
