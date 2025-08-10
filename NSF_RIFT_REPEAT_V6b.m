%%NSF_Rift_Repeat_V6
tic
clear all
close all

%%
%BASIC USER SETTINGS
Breath_Real=5000; % GIVEN Breath of domain in meters
S=1; %Breath Scaling-- 
Breath=S*Breath_Real; %Scaled Breath

cols=101; %Number of columns--equally spaced
Delx=Breath/(cols-1);
 
%Set Dely from Delx Through specifying an element Aspect Ratio (xdim/ydim)
AR=20;
Dely=Delx/AR;
%NOTE setting AR=n is essentially the same as setting S=1/n

sealevel=495; %sea-level elevation
consea=1; %salt concentration in sea 
rhorel=1.025;%relative density of saturated saline rhost/rhow

eps=0.35; %Porosity
storage=0.001; %Storage

totstep=5000; %5000; %XXXXXX5000; %Number of time steps of size 
delt=20; % Time step days
 
%hydraulic conductivity m/day
kxval=1;
kyval=.01;
%diffusion and dispersion coefficients
Dmol= 0.01;% m/day
aL=100;
aT=10; 

%XXXX USER DEFINES DOMAIN Here:
%By setting vector of x, y points on top surface 
xtop_user=[0,Breath];
ytop_user=[500,490];
% By setting vector of x, y points on bottom surface
ybot_user=[450,440];
xbot_user=[0,Breath];

%% SIMPLE MESH GENERATION--%USER INPUT NOT NEEDED

%We develop the full grid by interpolation in given user data
xq=linspace(0,Breath,cols);
ytop=interp1(xtop_user,ytop_user,xq);
ybot=interp1(xbot_user,ybot_user,xq);


%x and y values of nodes--counting up columns

kk=1;
Btop=[];
Btopsea=[];
for jj=1:cols
    %Try to fit rows on col to get as close as we can to target Dely
    rows=round(((ytop(jj)-ybot(jj))/Dely))+1;
    Delyreal=(ytop(jj)-ybot(jj))/(rows-1); %the actual spacing in a col.
    for ii=1:rows
        x(kk)=(jj-1)*Delx;
        y(kk)=ybot(jj)+(ii-1)*Delyreal;
        kk=kk+1;
    end
    ncol(jj)=rows; %number of nodes in column.
 
    %find Btop, Btopsea
    Btop=[Btop,kk-1];
    if ytop(jj)<sealevel
        Btopsea=[Btopsea,kk-1];
    end
    
end

N=size(x,2);

%Make mesh
t=[];
first=1;
for jj=1:cols-1
    last=first+ncol(jj)+ncol(jj+1)-1;
    xcc=x(first:last);
    
%     last - first + 1
    ycc=y(first:last);
    tccc=delaunay(xcc,ycc)
    tcc=delaunay(xcc,ycc)+first-1;
%     size(tcc)
    t=[t;tcc];
    first=first+ncol(jj);
end


Ntri=size(t,1);  %the size of t is the number of triangles

%max size of loose support
maxsup=0;
for ii=1:N
    maxsup=max(maxsup, size(find (t==ii),1));
end
maxsup=2*maxsup; %double the nodes in support

%printout grid
figure(1)
triplot(t,x,y, '-b');
%axis equal
hold on
pause(.1)
%END MESH XXXXXXXXXXXX
%% Geometric measures of grid--will need to be updated if grid chnages
%USER INPUT NOT REQUIRED
Volp=zeros(N,1);    %CV volume
xmid=zeros(Ntri,1); %elemnt mid point
ymid=zeros(Ntri,1);
Nx=zeros(Ntri,3);   %Derivatives of shape functions
Ny=zeros(Ntri,3);
tsup=zeros(N,1); %triangles in support
isup=ones(N,1);  %support index location
sup=ones(N,maxsup); %loose support


for itri=1:Ntri
    k1=t(itri,1); %global number of 1st node in trinagle itri
    k2=t(itri,2); %2nd node
    k3=t(itri,3); %3rd node
    v=(x(k2)*y(k3)-x(k3)*y(k2)-x(k1)*y(k3)+x(k1)*y(k2) ...
        +y(k1)*x(k3)-y(k1)*x(k2))/2;
    
    
    Volp(k1)=Volp(k1)+v/3; %contribution to control volume
    Volp(k2)=Volp(k2)+v/3;
    Volp(k3)=Volp(k3)+v/3;
    
    xmid(itri)=(x(k1)+x(k2)+x(k3))/3; %mid point of elemnt
    ymid(itri)=(y(k1)+y(k2)+y(k3))/3;
    
    %derivatives of shape functions
    Nx(itri,1)= (y(k2)-y(k3))/(2*v);   %Nx=shape function derivative
    Nx(itri,2)= (y(k3)-y(k1))/(2*v);   %the index 1, 2 or 3
    Nx(itri,3)= (y(k1)-y(k2))/(2*v);   %refers to local tri elemnt node
    Ny(itri,1)=-(x(k2)-x(k3))/(2*v);   %Ny=shape function derivative
    Ny(itri,2)=-(x(k3)-x(k1))/(2*v);
    Ny(itri,3)=-(x(k1)-x(k2))/(2*v);
   
    
    %loose support
    sup(k1,isup(k1))=k2;
    sup(k1,isup(k1)+1)=k3;
    isup(k1)=isup(k1)+2;
    tsup(k1)=tsup(k1)+1; 
    
    sup(k2,isup(k2))=k3;
    sup(k2,isup(k2)+1)=k1;
    isup(k2)=isup(k2)+2; 
    tsup(k2)=tsup(k2)+1; 
    
    sup(k3,isup(k3))=k1;
    sup(k3,isup(k3)+1)=k2;
    isup(k3)=isup(k3)+2;
    tsup(k3)=tsup(k3)+1; 
    
end

%% SET BOUNDARY CONDITIONS --USER INPUT NOT REQIIRED

%BOUNDARY NODE POINTS XXXXXXXXX
%Btop--stored node umbers on top boundary
%Btopsea--stored nodes at or below sea level.

%Arrays for Boundary conditions 
Big=1e18;
BCh=zeros(N,1); %boundary coeff values for head and solute
BCs=zeros(N,1);
BBh=zeros(N,1); %fixed head values
BBs=zeros(N,1); %fixed solute consentration values

%Fixed head value on top boundary 
%BC=Big, BB = value X Big
%top
BCh(Btop,1)=Big;
BBh(Btop,1)=Big*y(Btop);
BBh(Btopsea,1)=Big*(y(Btopsea)-(y(Btopsea)-sealevel)*(rhorel));  %correctinon for sea water 

BCs(Btop,1)=Big;
BBs(Btopsea)=consea*Big;

%%
%Transient Solution (implicit)--some user setting required
%Here we use an implicit solution-using a point iteration scheme
%To avoid re-calulation between iterations we assumme that, in head 
%calculation the value of the relative density is lagged one time step.
 
phi=max(ytop)*ones(N,1); %initial values maximum elevation
%phi=y; 
%con=zeros(N,1);
con=ones(N,1);
phinew=max(ytop)*ones(N,1);
%phinew=y;
%connew=zeros(N,1);
connew=ones(N,1);
qx=zeros(Ntri,1); %discharge
qy=zeros(Ntri,1);
con(Btopsea)=1;
connew(Btopsea)=1;


tolh=1e-3; %XXXXX LOOK AT THIS convergence tolarence on implcit iterations--head  XXXXHERE
tols=1e-7; %                                           --solute 
cntt = 0;
for tstep=1:totstep
%     tstep
    sto=storage;
    if tstep==1  %on firts step solve for steady flow XXXXXXX
        sto=0.0;
    end
    
    %Head Solution
    
    %coefficents
    ap=zeros(N,1);
    asup=zeros(N,maxsup);
    isup=ones(N,1);  %point in support 
   
   
    BBvar=zeros(N,1); %Variable density source
    
    
    kx=zeros(Ntri,1); %x direction conductivity value
    ky=zeros(Ntri,1); %y direction conductivity value
    
    for itri=1:Ntri

        %Code for elemnet conductivity 
           kx(itri)=kxval;
           ky(itri)=kyval; 
           
        
        
        cyc=[1,2,3;2,3,1;3,1,2];
        for node=1:3
            ii=cyc(node,1);
            jj=cyc(node,2);
            kk=cyc(node,3);
            
            k1=t(itri,ii); %global node number of element vertices
            k2=t(itri,jj);
            k3=t(itri,kk);
            
            Nx1=Nx(itri,ii);   %Nx=shape function derivative
            Nx2=Nx(itri,jj);
            Nx3=Nx(itri,kk);
            Ny1=Ny(itri,ii);   %Ny=shape function derivative
            Ny2=Ny(itri,jj);
            Ny3=Ny(itri,kk);
            
            
           %Face1
            delx= (x(k1)+x(k2)+x(k3))/3-(x(k1)+x(k2))/2;
            dely= (y(k1)+y(k2)+y(k3))/3-(y(k1)+y(k2))/2;
            
            face1_k1=S^2*kx(itri)*Nx1*dely-ky(itri)*Ny1*delx; %XXXX
            face1_k2=S^2*kx(itri)*Nx2*dely-ky(itri)*Ny2*delx;
            face1_k3=S^2*kx(itri)*Nx3*dely-ky(itri)*Ny3*delx;
      
            %varibale density source(face value)
            BBvar(k1)=BBvar(k1)-ky(itri)*((rhorel-1)/12)*(5*con(k1)+5*con(k2)+2*con(k3))*delx;
            
            
            %Face2
            delx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2;
            dely= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2;
         
            
            face2_k1=S^2*kx(itri)*Nx1*dely-ky(itri)*Ny1*delx; 
            face2_k2=S^2*kx(itri)*Nx2*dely-ky(itri)*Ny2*delx;
            face2_k3=S^2*kx(itri)*Nx3*dely-ky(itri)*Ny3*delx;

            
            ap(k1)=ap(k1)      -face1_k1-face2_k1;
            asup(k1,isup(k1))  =face1_k2+face2_k2; 
            asup(k1,isup(k1)+1)=face1_k3+face2_k3;
            isup(k1)=isup(k1)+2;
                       
            %varibale density source(face value)
            BBvar(k1)=BBvar(k1)-ky(itri)*((rhorel-1)/12)*(5*con(k1)+2*con(k2)+5*con(k3))*delx;
        end
    end
 
    %Update Head
  
    %implicit--used Jacobi (works well with vectorization)
    conver=1;
    phipre=phinew;
    
    iii = 0 ;
    while conver>tolh %convergence on iteration
       iii = iii + 1;
%        disp(iii);
       RHS=sum(asup.*phinew(sup),2);
       phinew=(sto*Volp.*phi+delt*(RHS+BBvar)+BBh)./(sto*Volp+BCh+ap*delt);
       conver=max(abs(phinew-phipre));  
       phipre=phinew;
    end
   
    
   
    phi=phinew; %old for new
    
%     fileID = fopen("aoutput.txt",'r');
%     formatSpec = '%f';
%     A = fscanf(fileID,formatSpec);
%     phi = A;

    %Solute solution
    ap=zeros(N,1);
    asup=zeros(N,maxsup);  %coeficents of lose unstructured support 
    isup=ones(N,1);   %number of nodes in lose unstructured support  
    
    
    for itri=1:Ntri
        
        %Diffusion
        cyc=[1,2,3;2,3,1;3,1,2];
        for node=1:3
            ii=cyc(node,1);
            jj=cyc(node,2);
            kk=cyc(node,3);
            
            k1=t(itri,ii); %global node number of element vertices
            k2=t(itri,jj);
            k3=t(itri,kk);
            
            Nx1=Nx(itri,ii);   %Nx=shape function derivative
            Nx2=Nx(itri,jj);
            Nx3=Nx(itri,kk);
            Ny1=Ny(itri,ii);   %Ny=shape function derivative
            Ny2=Ny(itri,jj);
            Ny3=Ny(itri,kk);
            
            
            if node==1 %discharge for element
                
 
              
                % --contribution from fresh-water head
                qxval=-S*kx(itri)*(Nx1*phi(k1)+Nx2*phi(k2)+Nx3*phi(k3));
                qyval=-ky(itri)*(Ny1*phi(k1)+Ny2*phi(k2)+Ny3*phi(k3));
                
                %store discharge for velcoity vector plot and D calc
                %these ARE in real sapce
                %Here they are caculated at element mid points
                %In setting up coeficicents (see later code)
                %they are calculated on face
                qx(itri)=qxval; 
                qy(itri)=qyval-ky(itri)*((rhorel-1)/3)*(con(k1)+con(k2)+con(k3));
 
                %Codeing for dispersion coeff --accounting for tensor
                %also asummes rel density contribution is constant in element
                 
                 qx2=qx(itri)^2;
                 qy2=qy(itri)^2;
                 qabs=sqrt(qx2+qy2);
                 
%                Dxx=(aL*qx2/qabs+aT*qy2/qabs+Dmol*eps);  
%                Dyy=aT*qx2/qabs+aL*qy2/qabs+Dmol*eps;  
%                Dxy=(aL-aT)*(qx(itri))*(qy(itri))/qabs;
                 
                 %Alt from Bear 1972 SAME AS above
                 Dxx=aT*qabs+(aL-aT)*qx2/qabs+Dmol*eps;
                 Dyy=aT*qabs+(aL-aT)*qy2/qabs+Dmol*eps;
                 Dxy=(aL-aT)*(qx(itri))*(qy(itri))/qabs;
      
             end
            
            %Face1
            qxface=S*qxval;
            qyface=qyval-ky(itri)*((rhorel-1)/12)*(5*con(k1)+5*con(k2)+2*con(k3));
            delx= (x(k1)+x(k2)+x(k3))/3-(x(k1)+x(k2))/2;
            dely= (y(k1)+y(k2)+y(k3))/3-(y(k1)+y(k2))/2; 
            
            face1_k1=(S^2*Dxx*Nx1+S*Dxy*Ny1)*dely-(Dyy*Ny1+S*Dxy*Nx1)*delx;
            face1_k2=(S^2*Dxx*Nx2+S*Dxy*Ny2)*dely-(Dyy*Ny2+S*Dxy*Nx2)*delx;
            face1_k3=(S^2*Dxx*Nx3+S*Dxy*Ny3)*dely-(Dyy*Ny3+S*Dxy*Nx3)*delx;
            
%             break
            %upwind
            qout=qxface*dely-qyface*delx; %flow out of vol k
            if qout>=0
                face1_k1=face1_k1-qout;
            else
                face1_k2=face1_k2-qout;
            end
            
            %Face2
            qxface=S*qxval;
            qyface=qyval-ky(itri)*((rhorel-1)/12)*(5*con(k1)+2*con(k2)+5*con(k3));
            
            delx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2;
            dely= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2;
            
            face2_k1=(S^2*Dxx*Nx1+S*Dxy*Ny1)*dely-(Dyy*Ny1+S*Dxy*Nx1)*delx;
            face2_k2=(S^2*Dxx*Nx2+S*Dxy*Ny2)*dely-(Dyy*Ny2+S*Dxy*Nx2)*delx;
            face2_k3=(S^2*Dxx*Nx3+S*Dxy*Ny3)*dely-(Dyy*Ny3+S*Dxy*Nx3)*delx;
      
            %upwind
            qout=qxface*dely-qyface*delx; %flow out of vol k
            if qout>=0
                face2_k1=face2_k1-qout;
            else
                face2_k3=face2_k3-qout;
            end
            
            %update diagonal element and support
            
            ap(k1)=ap(k1)      -face1_k1-face2_k1;
            asup(k1,isup(k1))  =face1_k2+face2_k2; 
            asup(k1,isup(k1)+1)=face1_k3+face2_k3;
            isup(k1)=isup(k1)+2;

         end
%          break
    end
    
%     break
%     disp(1234);
      
    %implicit--used Jacobi (works well with vectorization)
    conver=1;
    conpre=connew;
    cnt = 0;
    
    while conver>tols
        cnt = cnt + 1;
        RHS=sum(asup.*connew(sup),2);
        connew=(eps*Volp.*con+delt*RHS+BBs)./(eps*Volp+BCs+ap*delt);   
        conver=max(abs(connew-conpre));
        conpre=connew;
        
%         if cnt == 2;
%             break;
%         end
    end
    con=connew;
%     disp(max(con));
    
%     break
end

toc

quiver(xmid,ymid,qx,qy,'r','LineWidth',1)


figure (2)
trisurf(t,x,y,phi)
map= [0 0 1
    .3 .3 1
    .5 .5 1
 0 1 0
 0 .5 0
 .7 .7 0
 1 .4 0
 1 0 0];
colormap(map)
shading interp
view(0,90)
colorbar
%axis equal
hold on 


figure (3)
trisurf(t,x,y,con)
map= [0 0 1
    .3 .3 1
    .5 .5 1
 0 1 0
 0 .5 0
 .7 .7 0
 1 .4 0
 1 0 0];
 colormap(map)
shading interp
view(0,90)
colorbar
%axis equal
hold on 

dlmwrite('phi.txt',phi)

