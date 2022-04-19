function z = dranchuk(Pr,Tr,index)
%dranchuk subroutine to calcuate z values
%prepared by Yogesh Bansal
%this subroutine can work for calculating one single value of 'z' or for
%the entire matrix.
%if you want to calculate 'z' for one single block then you call this
%sub-routine as z = dranchuc(Pr,Tr,1) 
%for entire matrix, you need to send Pr and Tr for entire matrix with the
%index matrix.

A=[.3265; -1.07; -0.5339; 0.01569; -0.05165; 0.5475; -0.7361; 0.1844; ...
    0.1056; 0.6134; 0.721];

[u,v] = size(index);

c1 = A(1) + A(2)/Tr + A(3)*Tr^-3 + A(4)*Tr^-4 + A(5)*Tr^-5;
c2 = A(6) + A(7)/Tr + A(8)/Tr^2;
c3 = A(7)/Tr + A(8)/Tr^2;

for i=1:u
    for j=1:v
        if index(i,j) ~= 0
            ro(i,j) = 0.27*Pr(i,j)/Tr;
            diff = 1;
            while diff > 0.0001
                f = 1 - 0.27*Pr(i,j)/(ro(i,j)*Tr) + c1*ro(i,j) + ...
                    c2*ro(i,j)^2 - A(9)*c3*ro(i,j)^5 + ... 
                    A(10)*(1 + A(11)*ro(i,j)^2)*ro(i,j)^2*exp(-1*A(11)*ro(i,j)^2)/Tr^3;
                fd = 0.27*Pr(i,j)/(ro(i,j)^2*Tr) + c1 + 2*c2*ro(i,j) - ...
                    5*A(9)*c3*ro(i,j)^4 + ... 
                    2*A(10)*ro(i,j)*(1+A(11)*ro(i,j)^2 - ...
                    A(11)^2*ro(i,j)^4)*exp(-1*A(11)*ro(i,j)^2)/Tr^3;
                
                ron(i,j) = ro(i,j) - f/fd;
                diff = abs(ron(i,j) - ro(i,j));
                ro(i,j) = ron(i,j);
            end
            z(i,j) = 0.27*Pr(i,j)/(ro(i,j)*Tr);
        end
    end
end


