%%NSF_Delta_pro_production_V1
%Simulated advanceing simplified clinoform-Field Scale
%assumed initial condition
%revised way of tacking salt trapping
%Bangledesh settings

clear all
close all

%%
%BASIC USER SETTINGS

%Basic Geometric data
Breath=300000; %%Breath (m) of domain (Length of long profile) 
cols=101; %Number of grid columns--equally spaced
Delx=Breath/(cols-1); %fixed width between columns
xcols=linspace(0,Breath,cols); %xlocations of columns
%Set Dely (Typical height of an elemnet) -
% from Delx, through specifying an element aspect Ratio (xdim/ydim)
AR=80;
Dely=Delx/AR;

delt=100; % Time step days

%Physical properties
consea=1; %Salt concentration in sea 
rhorel=1.025;%Relative density of saturated saline rhosat/rhowater

eps=0.35; %Aquifer porosity
sto=0.001; %Aquifer storage


%Hydraulic conductivity m/day
kxval=10;
kyval=kxval/100;
%diffusion and dispersion coefficients
Dmol= 0.00001;% m^2/day
aL=50;  % m keep here for now
aT=aL/10; 

%% USER Defined synthetic stratigraphy
%models clinerform as assumming consrant fluivila and foreset slopes
sealevel=300; %600; %Sea-level setting
basement=250; % 500 Depth of basement
foreslope=.05; %foreset slope 05
topslope=0.0005; %0.0005; %fluival topset slope

%initila position (m) of thelocation of toe 
%(foreset interection with basemnet) 
xtoein=(sealevel-basement)/foreslope; %Assuming NO intial fluvial surface
xtoe=xtoein;

toevel=.04;  %toe velocity over basement m/day

%Creatate and store evloving margin and clinoform
tstep=1;
while xtoe<0.8*Breath
    xtoe=xtoein+toevel*(tstep-1)*delt;
    xshore=xtoe-(sealevel-basement)/foreslope;
    eta(tstep,:)=min ( max((xshore-xcols)*topslope+sealevel,sealevel), ...
                       max((xtoe-xcols)*foreslope+basement,basement) );
    etabot(tstep,:)=zeros(1,cols);
    tstep=tstep+1;
end
totstep=tstep-1; %total timesteps of calculation-time for tor to get to 90%


%INITIAL DOMAIN 
ytop=eta(1,:);
ybot=etabot(1,:);
ncol=round((ytop-ybot)/Dely)+1;
N=sum(round((ytop-ybot)/Dely)+1);
phinew=max(ytop)*ones(N,1);

%set salt intitial cons of points vertically below sea level to value 1.
connew=[];
for jj=1:cols
    if ytop(jj)>sealevel
        connew=[connew;zeros(ncol(jj),1)];
    else
        connew=[connew;ones(ncol(jj),1)];
    end
end
       

%trapping number - small=not much trapping
TR=(sealevel/aL)*(eps*toevel/kxval/topslope) %(sealevel/aL) is peclet-lower val of peclet number is more dispersion? controlling in early time but becomes steady state in later time and other half of the eq dominates

%% Time Stepping

tolh=1e-7; %convergence tolerance
tols=1e-7;                                      

printval=50; 
printspace=50; % time slice interval
ipr = 0; % time slice for tecplot output

% dimension plot variables for tecplot output
%  MAP

NPR = totstep/printspace; 
NNPR = round(NPR);
MaxNnode = 2000;
MaxNTri = 4000;
xarray = zeros(6,MaxNnode,NNPR);
elarray = zeros (3,MaxNTri,NNPR);

totstep = 1
aaa = zeros(totstep,20)
for tstep=1:totstep
    tstep
    
%% Adjustmnent of strat and set sea-level

ytop=eta(tstep,:);
ybot=etabot(tstep,:);


%x and y values of nodes--counting up columns
ncolold=ncol;
ncol=round((ytop-ybot)/Dely)+1;

kk=1;
Btop=[]; %list of nodes on top surface
Btopsea=[];  %list of top domain nodes under sea  
x=[];
y=[];
for jj=1:cols
    x=[x;Delx*(jj-1)*ones(ncol(jj),1)];
    y=[y;linspace(ybot(jj),ybot(jj)+(ncol(jj)-2)*Dely,ncol(jj)-1)';ytop(jj)];
    %Find Btop, Btopsea Nodes on top surface, topsurface nodes under sea
    Btop=[Btop,size(y,1)];
    if ytop(jj)<sealevel
        Btopsea=[Btopsea,size(y,1)];
    end
end
Bbot=[1,Btop(1:cols-1)+1]; %List of bottom nodes

N=size(x,1);
%resize phi and con as grid expands 
con=[];
phi=[];
last=0;
for jj=1:cols
    first=last+1;
    last=first+ncolold(jj)-1;
    
    if ncolold(jj)==ncol(jj)
        phi=[phi;phinew(first:last)];
        con=[con;connew(first:last)];
    end
    
    if ncolold(jj)>ncol(jj)  %errosion
        phi=[phi;phinew(first:last-1)];
        con=[con;connew(first:last-1)];
    end
    
    if ncolold(jj)<ncol(jj)  %deposition
        phi=[phi;phinew(first:last);phinew(last)];
        con=[con;connew(first:last);connew(last)];
        %interpolate to find value at added point
        mnode=size(con,1);
%         con(mnode)
%         con(mnode-1)
%         phi(mnode)
%         phi(mnode-1)
        yrat=(y(mnode)-y(mnode-1))/(y(mnode)-y(mnode-2));
%         phi(mnode)
        phi(mnode-1)=phi(mnode)-yrat*(phi(mnode)-phi(mnode-2));
        con(mnode-1)=con(mnode)-yrat*(con(mnode)-con(mnode-2));
%         break
    end   
end


phinew=phi;
connew=con;

%% Make Mesh

t=[];
first=1;
for jj=1:cols-1
    last=first+ncol(jj)+ncol(jj+1)-1;
    xcc=x(first:last);
    ycc=y(first:last);
    tcc=delaunay(xcc,ycc)+first-1;
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
% maxsup = 20

%printout initial grid 
if tstep==1
figure(1)
triplot(t,x,y, '-b');
%hold on
%plot(polyshape([50000 50000 250000 250000],[175 100 100 175])); %lowcondunit
%axis equal
end


%END MESH XXXXXXXXXXXX

%% Geometric properties of mesh 
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


%in_mat1=inpolygon(xmid,ymid,[50000 50000 250000 250000],[175 100 100 175]); %might need to move to loop bc strat changes


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
%correctinon for sea water 
BBh(Btopsea,1)=Big*(y(Btopsea)-(y(Btopsea)-sealevel)*(rhorel));  
BCs(Btop,1)=Big;
BBs(Btopsea)=consea*Big;

%% Head Solution

    %coefficents
    ap=zeros(N,1);
    asup=zeros(N,maxsup);
    isup=ones(N,1);  %point in support 
   
    BBvar=zeros(N,1); %Variable density source
    
    kx=zeros(Ntri,1); %x direction conductivity value
    ky=zeros(Ntri,1); %y direction conductivity value
    
    for itri=1:Ntri

        %Code for elemnet conductivity 
        %if(in_mat1(itri)==1)
            %kx(itri)=kxval/1000;
            %ky(itri)=kyval/1000;
        %else
            kx(itri)=kxval;
            ky(itri)=kyval; 
        %end
  
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
            
            face1_k1=kx(itri)*Nx1*dely-ky(itri)*Ny1*delx;
            face1_k2=kx(itri)*Nx2*dely-ky(itri)*Ny2*delx;
            face1_k3=kx(itri)*Nx3*dely-ky(itri)*Ny3*delx;
      
            %varibale density source(face value)
            BBvar(k1)=BBvar(k1)-ky(itri)*((rhorel-1)/12)*(5*con(k1)+5*con(k2)+2*con(k3))*delx;
            
            
            %Face2
            delx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2;
            dely= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2;
         
            
            face2_k1=kx(itri)*Nx1*dely-ky(itri)*Ny1*delx; 
            face2_k2=kx(itri)*Nx2*dely-ky(itri)*Ny2*delx;
            face2_k3=kx(itri)*Nx3*dely-ky(itri)*Ny3*delx;

            
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
    
    while conver>tolh %convergence on iteration
       RHS=sum(asup.*phinew(sup),2);
       phinew=(sto*Volp.*phi+delt*(RHS+BBvar)+BBh)./(sto*Volp+BCh+ap*delt);
       conver=max(abs(phinew-phipre));  
       phipre=phinew;
    end
   
    phi=phinew; %old for new
    
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
            
            
            if node==1 %Dispersion for element
              
                % --contribution to discharge from fresh-water head
                qxval=-kx(itri)*(Nx1*phi(k1)+Nx2*phi(k2)+Nx3*phi(k3));
                qyval=-ky(itri)*(Ny1*phi(k1)+Ny2*phi(k2)+Ny3*phi(k3));
                
                %actual discharges at elemnt midpoints 
                %used to calculate dispersion tensor
                %also asummes rel density contribution is constant in element
 
                qxmid=qxval; 
                qymid=qyval-ky(itri)*((rhorel-1)/3)*(con(k1)+con(k2)+con(k3));
  
            qx(itri) = qxmid;       % MAP
            qy(itri) = qymid;       % MAP
            
                qx2=qxmid^2;
                qy2=qymid^2;
                qabs=sqrt(qx2+qy2);    
                 
                %Dispersion from Bear 1972 SAME AS above
                Dxx=aT*qabs+(aL-aT)*qx2/qabs+Dmol*eps;
                Dyy=aT*qabs+(aL-aT)*qy2/qabs+Dmol*eps;
                Dxy=(aL-aT)*(qxmid)*(qymid)/qabs;
      
             end
             
            %Face1
            qxface=qxval;
            qyface=qyval-ky(itri)*((rhorel-1)/12)*(5*con(k1)+5*con(k2)+2*con(k3));
            delx= (x(k1)+x(k2)+x(k3))/3-(x(k1)+x(k2))/2;
            dely= (y(k1)+y(k2)+y(k3))/3-(y(k1)+y(k2))/2; 
            
            face1_k1=(Dxx*Nx1+Dxy*Ny1)*dely-(Dyy*Ny1+Dxy*Nx1)*delx;
            face1_k2=(Dxx*Nx2+Dxy*Ny2)*dely-(Dyy*Ny2+Dxy*Nx2)*delx;
            face1_k3=(Dxx*Nx3+Dxy*Ny3)*dely-(Dyy*Ny3+Dxy*Nx3)*delx;
            
   
            %upwind
            qout=qxface*dely-qyface*delx; %flow out of vol k
            if qout>=0
                face1_k1=face1_k1-qout;
            else
                face1_k2=face1_k2-qout;
            end
            
            %Face2
            qxface=qxval;
            qyface=qyval-ky(itri)*((rhorel-1)/12)*(5*con(k1)+2*con(k2)+5*con(k3));

            
            delx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2;
            dely= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2;
            
            face2_k1=(Dxx*Nx1+Dxy*Ny1)*dely-(Dyy*Ny1+Dxy*Nx1)*delx;
            face2_k2=(Dxx*Nx2+Dxy*Ny2)*dely-(Dyy*Ny2+Dxy*Nx2)*delx;
            face2_k3=(Dxx*Nx3+Dxy*Ny3)*dely-(Dyy*Ny3+Dxy*Nx3)*delx;
      
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
         
      end
    
      
    %implicit--used Jacobi (works well with vectorization)
    conver=1;
    conpre=connew;
    while conver>tols
        RHS=sum(asup.*connew(sup),2);
        connew=(eps*Volp.*con+delt*RHS+BBs)./(eps*Volp+BCs+ap*delt);   
        conver=max(abs(connew-conpre));
        conpre=connew;
    end
    
 
    %connew=min(connew,1); %XXXX mass Limiter
    con=connew;
    
    %% AGE solution

     %age solution
     age = zeros(809,1);
     agenew = age;

     Q=zeros(N,1);

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
             
             
             if node==1 %Dispersion for element
               
                 % --contribution to discharge from fresh-water head
                 qxval=-kx(itri)*(Nx1*phi(k1)+Nx2*phi(k2)+Nx3*phi(k3));
                 qyval=-ky(itri)*(Ny1*phi(k1)+Ny2*phi(k2)+Ny3*phi(k3));
                 
                 %actual discharges at elemnt midpoints 
                 %used to calculate dispersion tensor
                 %also asummes rel density contribution is constant in element
  
                 qxmid=qxval; 
                 qymid=qyval-ky(itri)*((rhorel-1)/3)*(age(k1)+age(k2)+age(k3));
   
             qx(itri) = qxmid;       % MAP
             qy(itri) = qymid;       % MAP
             
                 qx2=qxmid^2;
                 qy2=qymid^2;
                 qabs=sqrt(qx2+qy2);    
                  
                 %Dispersion from Bear 1972 SAME AS above
                 Dxx=aT*qabs+(aL-aT)*qx2/qabs+Dmol*eps;
                 Dyy=aT*qabs+(aL-aT)*qy2/qabs+Dmol*eps;
                 Dxy=(aL-aT)*(qxmid)*(qymid)/qabs;
       
              end
              
             %Face1
             qxface=qxval;
             qyface=qyval-ky(itri)*((rhorel-1)/12)*(5*age(k1)+5*age(k2)+2*age(k3));
             delx= (x(k1)+x(k2)+x(k3))/3-(x(k1)+x(k2))/2;
             dely= (y(k1)+y(k2)+y(k3))/3-(y(k1)+y(k2))/2; 
             
             face1_k1=(Dxx*Nx1+Dxy*Ny1)*dely-(Dyy*Ny1+Dxy*Nx1)*delx;
             face1_k2=(Dxx*Nx2+Dxy*Ny2)*dely-(Dyy*Ny2+Dxy*Nx2)*delx;
             face1_k3=(Dxx*Nx3+Dxy*Ny3)*dely-(Dyy*Ny3+Dxy*Nx3)*delx;

             Q(k1)=Q(k1)-Volp(k1);
             
    
             %upwind
             qout=qxface*dely-qyface*delx; %flow out of vol k
             if qout>=0
                 face1_k1=face1_k1-qout;
             else
                 face1_k2=face1_k2-qout;
             end
             
             %Face2
             qxface=qxval;
             qyface=qyval-ky(itri)*((rhorel-1)/12)*(5*age(k1)+2*age(k2)+5*age(k3));
 
             
             delx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2;
             dely= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2;
             
             face2_k1=(Dxx*Nx1+Dxy*Ny1)*dely-(Dyy*Ny1+Dxy*Nx1)*delx;
             face2_k2=(Dxx*Nx2+Dxy*Ny2)*dely-(Dyy*Ny2+Dxy*Nx2)*delx;
             face2_k3=(Dxx*Nx3+Dxy*Ny3)*dely-(Dyy*Ny3+Dxy*Nx3)*delx;
       
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
          
       end
     
       
     %implicit--used Jacobi (works well with vectorization)
     conver=1;
     agepre=agenew;
     while conver>tols
         RHS=sum(asup.*agenew(sup),2);
         agenew=(eps*Volp.*age+delt*(RHS+Q)+BBs)./(eps*Volp+BCs+ap*delt);   
         conver=max(abs(agenew-agepre));
         agepre=agenew;
     end
     
  
     %connew=min(connew,1); %XXXX mass Limiter
     age=agenew;
age
 %% END of AGE solution

%print out at set time steps    
if tstep>printval
    ipr = ipr+1;
    printval=printval+printspace;
    figure (2)
    plot(xcols,ytop,'k-', 'LineWidth',2)
    %title('water less that 10% saline')
    ylim([0 450])
    hold on
    trisurf(t,x,y,con)
 
 %%%%%%%%%% Store Tecplot variables for later plotting %%%%%%  MAP       
        for itri=1:Ntri
            vxp(t(itri,1))= qx(itri);
            vxp(t(itri,2))= qx(itri);
            vxp(t(itri,3))= qx(itri);
            
            vyp(t(itri,1))= qy(itri);
            vyp(t(itri,2))= qy(itri);
            vyp(t(itri,3))= qy(itri);   
        end
        nelem(ipr) = Ntri;
        nnode(ipr) = N;
            
    for n=1:N
        xarray(1,n,ipr)=x(n);
        xarray(2,n,ipr)=y(n);
        xarray(3,n,ipr)=phi(n);
        xarray(4,n,ipr)=con(n);
        xarray(5,n,ipr)=vxp(n);
        xarray(6,n,ipr)=vyp(n);        
    end

    for itri=1:Ntri
        if (t(itri,1)>0)
        elarray(1,itri,ipr)=t(itri,1);
        elarray(2,itri,ipr)=t(itri,2);
        elarray(3,itri,ipr)=t(itri,3);
        end
    end
    
    
    
    map=[0 0 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1];
    colormap(map)
    shading interp
    view(0,90)
    daspect([25 1 1])
    txt= ['$K_x =$' num2str(kxval) ', $v_{toe} =$' num2str(toevel) ...
          ', topset =' num2str(topslope) ', foreset ='  num2str(foreslope) ...
             ', TR =' num2str(TR) ];
    text(5000,390,txt,'fontsize', 14,'interpreter','latex')
    
    hold off
    
    % find con .5 location on base--wait till corner point is C>.4
    if con(1)<0.4
        conbotvals=con(Bbot);
        xq=x(Bbot);
        inval=interp1(conbotvals-(1e-10)*xq,xq,.5);
        tval=delt*tstep;
        shpos=tval*toevel;
        toe_shor(ipr) = shpos;
        toe_salt(ipr) = inval;
        
        %difference between shore and con=0.5
        shpos-inval %lag
        figure (11)
        title('movment of shoreline and C=0.5 isoline')
        plot(tval,shpos,'k*')
        hold on
        plot(tval,inval,'k.')
        legend('shoreline', 'C=0.5 isoline')
        
        
    end
    pause(.1)
    
end
    
    
    
end





figure (3)
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
daspect([25 1 1])



figure (4)
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
%daspect([25 1 1])



   figure (5)
   subplot(2,1,2)
    plot(xcols,ytop,'k-', 'LineWidth',2)
    %title('C, TR=0.1')
    ylim([0 450])
    hold on
    trisurf(t,x,y,con)
    
    map=[0 0 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1
        1 1 1];
    colormap(map)
    shading interp
    view(0,90)
    
    txt= ['$K_x =$' num2str(kxval) ', $v_{toe} =$' num2str(toevel)...
          ', topset =' num2str(topslope) ', foreset ='  num2str(foreslope)];
    text(5000,390,txt,'fontsize', 14,'interpreter','latex')
    hold on 
   
figure(10)  
%%%%%%%%%%%%%%%%% Tecplot Output %%%%%%%%%%%%%%%%%%%  MAP
output_tec_filename='voller_tec.dat';
fid=fopen(output_tec_filename,'w');
fprintf(fid,'variables = "x", "z","head","conc","vx","vz"\n');



for ipr=1:NPR
    
       nn = nnode(ipr);
       ne = nelem(ipr);
        xarray2 = xarray(:,1:nn,ipr);
        elarray2 = elarray(:,1:ne,ipr);


fprintf(fid,'ZONE N= %d, E= %d, F=FEPOINT, ET=TRIANGLE\n',nn,ne);
fprintf(fid,'%6.3f %6.3f %6.3f %6.3f %15.8f %15.8f\n',xarray2);
fprintf(fid,'%d %d %d\n',elarray2);
end
fclose(fid);

figure(6)
quiver(xmid,5*ymid,qx',qy','-k','MaxHeadSize',0.001)
