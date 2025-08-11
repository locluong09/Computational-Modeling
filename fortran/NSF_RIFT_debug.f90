program nsf_rift
	implicit none

	real :: breath_real = 5000.0 ! given breath domain in meters
	real :: s = 1.0 ! breath scaling
	real :: breath ! scaled Breath

	real :: delx, dely
	real :: ar = 20.0
	real :: sealevel = 495.0 ! sea level evaluation
	real :: consea = 1.0 ! salt concentration in sea
	real :: rhorel = 1.025 ! relative density of saturated saline
	real :: eps = 0.35 ! porosity
	real :: storage = 0.001
	real :: sto

	integer :: cols = 101, rows = 21 ! number of columns and rows
	integer :: totstep = 5000 ! number of time steps
	integer :: delt = 20
	integer :: tstep

	! hydraulic conductivity m/day
	real :: kxval = 1.0
	real :: kyval = 0.01
	! diffusion and dispersion coefficients
	real :: Dmol = 0.01
	real :: aL = 100.0
	real :: aT = 10.0

	integer :: idx
	integer :: i =1, j=1, k =1
	integer :: ii, jj, kk
	integer :: maxsup = 12
	integer :: N, Ntri
	integer :: iii = 0 !support Btopsea

	! USER SETTING 
	real, dimension(1:2) :: xtop_user !vector of x, y points on top surface 
	real, dimension(1:2) :: ytop_user 
	real, dimension(1:2) :: xbot_user !vector of x, y points on bottom surface
	real, dimension(1:2) :: ybot_user 

	integer, dimension(1:101)	 :: xq
	real, dimension(1:101)	 :: ytop
	real, dimension(1:101)	 :: ybot

	real :: Delyreal

	integer, dimension(1:101) :: Btop
	integer, dimension(1:50) :: Btopsea
	real, dimension(1:2121) :: x, y
	real, dimension(1:101) :: ncol
	
	integer, dimension(3,4000) :: temp
	integer, dimension(4000,3) :: t

	integer :: k1,k2,k3

	real, allocatable, dimension(:) :: Volp
	real, allocatable, dimension(:) :: xmid
	real, allocatable, dimension(:) :: ymid
	real, allocatable, dimension(:,:) :: Nx
	real, allocatable, dimension(:,:) :: Ny

	real, allocatable, dimension(:) :: BCh
	real, allocatable, dimension(:) :: BCs
	real, allocatable, dimension(:) :: BBh
	real, allocatable, dimension(:) :: BBs

	integer, allocatable, dimension(:) :: tsup ! triangle in support
	integer, allocatable, dimension(:) :: isup ! support index location
	integer, allocatable, dimension(:,:) :: sup !loose support
	real, allocatable, dimension(:,:) :: temp_sup ! use for temporary multiplication

	real :: v
	real :: Big = 1e18

	real, allocatable, dimension(:) :: phi
	real, allocatable, dimension(:) :: con
	real, allocatable, dimension(:) :: phinew
	real, allocatable, dimension(:) :: phipre
	real, allocatable, dimension(:) :: connew
	real, allocatable, dimension(:) :: conpre
	real, allocatable, dimension(:) :: qx
	real, allocatable, dimension(:) :: qy

	real :: tolh = 1e-7 ! convergence tolarence for head
	real :: tols = 1e-7 ! convergence tolarence for solute

	real, allocatable, dimension(:) :: ap
	real, allocatable, dimension(:,:) :: asup
	real, allocatable, dimension(:) :: BBvar
	real, allocatable, dimension(:) :: kx
	real, allocatable, dimension(:) :: ky

	integer, dimension(3,3) :: cyc
	integer :: node
	integer :: itri

	real :: Nx1, Nx2,  Nx3
	real :: Ny1, Ny2,  Ny3
	real :: dx, dy
	real :: face1_k1, face1_k2, face1_k3
	real :: face2_k1, face2_k2, face2_k3
	real :: conver = 1.0
	real :: qxval, qyval, qxface, qyface, qabs, qout
	real :: qx2, qy2
	real :: Dxx, Dyy, Dxy

	real, allocatable, dimension(:) :: RHS
	integer :: cnt = 0



	breath = s*breath_Real
	delx = breath/(cols-1)
	dely = delx/ar
	xtop_user = (/0.0, breath/)
	ytop_user = (/500.0, 490.0/)
	ybot_user = (/450.0, 440.0/)
	xbot_user = (/0.0, breath/)

	xq = (/(idx, idx = 0, 5000, 50)/)

	! fully grid by interpolating data
	ytop = interpolate1(xtop_user,ytop_user,xq)
	ybot = interpolate1(xbot_user,ybot_user,xq)

	do j = 1, cols
		Delyreal = (ytop(j)-ybot(j))/(rows-1)
		do i = 1, rows
			x(k) = (j-1)*Delx
			y(k)=ybot(j)+(i-1)*Delyreal
			k = k + 1
		end do
		ncol(j) = rows
		! find Btop, Btopsea
		Btop(j) = k - 1
		if (ytop(j) < sealevel) then
			iii = iii + 1
			Btopsea(iii) = k-1
		end if
	end do


	! t.txt contains mesh grid of trianglation
	open(unit=11, file="t.txt") 
	read(11, *) temp
	close(11)
	t = transpose(temp)


	Ntri = size(t,1) !the size of t is the number of triangles
	N = size(x)

	! allocate varibales

	allocate(Volp(N)) 
	allocate(xmid(Ntri)) 
	allocate(ymid(Ntri)) 
	allocate(Nx(Ntri,3))
	allocate(Ny(Ntri,3))
	allocate(tsup(N))
	allocate(isup(N))
	allocate(sup(N, maxsup))
	allocate(temp_sup(N, maxsup))

	allocate(BCh(N)) 
	allocate(BCs(N)) 
	allocate(BBh(N)) 
	allocate(BBs(N)) 

	allocate(phi(N))
	allocate(con(N))
	allocate(phinew(N))
	allocate(connew(N))
	allocate(qx(Ntri))
	allocate(qy(Ntri))

	allocate(ap(N))
	allocate(asup(N, maxsup))
	allocate(BBvar(N))
	allocate(kx(Ntri))
	allocate(ky(Ntri))
	allocate(phipre(N))
	allocate(conpre(N))

	allocate(RHS(N))

	! set zeros and ones arrays
	do i = 1,N
		Volp(i) = 0.0
		BCh(i) = 0.0
		BCs(i) = 0.0
		BBh(i) = 0.0
		BBs(i) = 0.0
		phi(i) = 1.0
		con(i) = 1.0
		connew(i) = 1.0
		phinew(i) = 1.0
		phipre(i) = 1.0
		BBvar(i) = 0.0
		RHS(i) = 0.0
		tsup(i) = 0
		isup(i) = 1
	end do

	do i = 1,Ntri
		xmid(i) = 0.0
		ymid(i) = 0.0
		qx(i) = 0.0
		qy(i) = 0.0
		kx(i) = 0.0
		ky(i) = 0.0
	end do

	do i = 1,Ntri
		do j = 1,3
			Nx(i,j) = 0.0
			Ny(i,j) = 0.0
		end do
	end do

	do i = 1,N
		do j = 1,maxsup
			sup(i,j) = 1
			asup(i,j) = 0.0
		end do
	end do


	! Find the derivative of shape function
	! and neighbor nodes around a support node
	do i = 1, Ntri
		k1 = t(i,1)
		k2 = t(i,2)
		k3 = t(i,3)
		v=(x(k2)*y(k3)-x(k3)*y(k2)-x(k1)*y(k3)+x(k1)*y(k2) &
        +y(k1)*x(k3)-y(k1)*x(k2))/2

        Volp(k1) = Volp(k1) + v/3
        Volp(k2) = Volp(k2) + v/3
        Volp(k3) = Volp(k3) + v/3

        xmid(i)=(x(k1)+x(k2)+x(k3))/3
    	ymid(i)=(y(k1)+y(k2)+y(k3))/3

    	! derivatives of shape functions
    	Nx(i,1)= (y(k2)-y(k3))/(2*v)
    	Nx(i,2)= (y(k3)-y(k1))/(2*v)
    	Nx(i,3)= (y(k1)-y(k2))/(2*v)
    	Ny(i,1)=-(x(k2)-x(k3))/(2*v)
    	Ny(i,2)=-(x(k3)-x(k1))/(2*v)
    	Ny(i,3)=-(x(k1)-x(k2))/(2*v)

    	! loose support
    	sup(k1,isup(k1))=k2
    	sup(k1,isup(k1)+1)=k3
    	isup(k1)=isup(k1)+2
    	tsup(k1)=tsup(k1)+1 
    
    	sup(k2,isup(k2))=k3
    	sup(k2,isup(k2)+1)=k1
    	isup(k2)=isup(k2)+2 
    	tsup(k2)=tsup(k2)+1 
    
    	sup(k3,isup(k3))=k1
    	sup(k3,isup(k3)+1)=k2
    	isup(k3)=isup(k3)+2
    	tsup(k3)=tsup(k3)+1
    end do 

    ! Set boundary conditions
    ! Btop--stored node umbers on top boundary
    ! Btopsea--stored nodes at or below sea level
    ! Fixed head value on top boundary 
	! BC=Big, BB = value X Big
    do i = 1, size(Btop)
    	j = Btop(i)
    	BCh(j) = Big
    	BBh(j) = Big*y(j)
    	BCs(j)=Big
    end do


    do i = 1, size(Btopsea)
    	j = Btopsea(i)
    	BBh(j)=Big*(y(j)-(y(j)-sealevel)*(rhorel))
    	BBs(j)=Big*consea
    end do

    phi = maxval(ytop)*phi
    phinew = maxval(ytop)*phinew

    do i = 1,size(Btopsea)
    	j = Btopsea(i)
    	con(j) = 1.0
    	connew(j) = 1.0
    end do


    
    ! Run simulation
    do tstep = 1, totstep

    	cnt = cnt + 1
    	! print out tstep for each iteration
    	write (6,*) cnt

    	sto = storage
    	if (tstep == 1) then
    		sto = 0.0
    	end if

    	! Head solution

    	do i = 1, N
    		ap(i) = 0.0
    		do j = 1, maxsup
    			asup(i,j) = 0.0
    		end do
    		BBvar(i) = 0.0
    		isup(i) = 1
    	end do

    	do i = 1, Ntri
    		kx(i) = 0.0
    		ky(i) = 0.0
    	end do

    	do itri = 1, Ntri
    		
    		kx(itri) = kxval
    		ky(itri) = kyval

    		! if (xmid(itri)>1500.0 && xmid(itri)<3500.0) .and. (ymid(itri)>460.0 && ymid(itri)<480.0) then
    		! 	kx(itri)=0.000001
    		! 	ky(itri)=0.000001
    		! end if

    		cyc = reshape((/ 1, 2, 3, 2, 3, 1, 3, 1, 2 /), shape(cyc))

    		do node = 1,3
    			ii=cyc(node,1)
            	jj=cyc(node,2)
            	kk=cyc(node,3)

            	
            	

            	k1=t(itri,ii)
            	k2=t(itri,jj)
            	k3=t(itri,kk)
            	
            	

            	

            	Nx1=Nx(itri,ii)
            	Nx2=Nx(itri,jj)
            	Nx3=Nx(itri,kk)
            	Ny1=Ny(itri,ii)
            	Ny2=Ny(itri,jj)
            	Ny3=Ny(itri,kk)

            	! Face1
            	dx= (x(k1)+x(k2)+x(k3))/3-(x(k1)+x(k2))/2
            	dy= (y(k1)+y(k2)+y(k3))/3-(y(k1)+y(k2))/2
            	

            	face1_k1=S**2*kx(itri)*Nx1*dy-ky(itri)*Ny1*dx
            	face1_k2=S**2*kx(itri)*Nx2*dy-ky(itri)*Ny2*dx
            	face1_k3=S**2*kx(itri)*Nx3*dy-ky(itri)*Ny3*dx

            	BBvar(k1)=BBvar(k1)-ky(itri)*((rhorel-1)/12)*(5*con(k1)+5*con(k2)+2*con(k3))*dx

            	! Face2
            	dx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2
            	dy= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2
            	
         
            	! print *, dx

            	face2_k1=S**2*kx(itri)*Nx1*dy-ky(itri)*Ny1*dx
            	face2_k2=S**2*kx(itri)*Nx2*dy-ky(itri)*Ny2*dx
            	! print *, face2_k1
            	face2_k3=S**2*kx(itri)*Nx3*dy-ky(itri)*Ny3*dx


            	ap(k1)= ap(k1) - face1_k1 - face2_k1
            	asup(k1,isup(k1))   = face1_k2 + face2_k2
            	asup(k1,isup(k1)+1) = face1_k3 + face2_k3
            	isup(k1)=isup(k1)+2

            	BBvar(k1)=BBvar(k1)-ky(itri)*((rhorel-1)/12)*(5*con(k1)+2*con(k2)+5*con(k3))*dx
            end do
    	end do

    	
    	! Update head
    	conver=1.0
    	phipre=phinew
    	do while (conver > tolh)
    		do i = 1,N
    			do j = 1,maxsup
    				temp_sup(i,j) = phinew(sup(i,j))
    			end do
    		end do
    		RHS    = sum(asup*temp_sup,2)
       		phinew = (sto*(Volp*phi)+delt*(RHS+BBvar)+BBh)/(sto*Volp+BCh+ap*delt)
       		conver = maxval(abs(phinew-phipre))
       		phipre = phinew


    	end do

    	phi=phinew
    	
    	
    		
    	

    	! Solute solution
    	do i = 1, N
    		ap(i) = 0.0
    		do j = 1, maxsup
    			asup(i,j) = 0.0
    		end do
    		isup(i) = 1
    	end do

    	do itri = 1, Ntri
    		cyc = reshape((/ 1, 2, 3, 2, 3, 1, 3, 1, 2 /), shape(cyc))
    		do node = 1,3
    			ii=cyc(node,1)
            	jj=cyc(node,2)
            	kk=cyc(node,3)

            	k1=t(itri,ii)
            	k2=t(itri,jj)
            	k3=t(itri,kk)

            	Nx1=Nx(itri,ii)
            	Nx2=Nx(itri,jj)
            	Nx3=Nx(itri,kk)
            	Ny1=Ny(itri,ii)
            	Ny2=Ny(itri,jj)
            	Ny3=Ny(itri,kk)

            	if (node == 1) then
            		qxval=-S*kx(itri)*(Nx1*phi(k1)+Nx2*phi(k2)+Nx3*phi(k3))
                	qyval=-ky(itri)*(Ny1*phi(k1)+Ny2*phi(k2)+Ny3*phi(k3))

                	qx(itri)=qxval
                	qy(itri)=qyval-ky(itri)*((rhorel-1)/3)*(con(k1)+con(k2)+con(k3))

                	qx2=qx(itri)**2
                 	qy2=qy(itri)**2
                 	qabs=sqrt(qx2+qy2)

                 	Dxx=aT*qabs+(aL-aT)*qx2/qabs+Dmol*eps
                 	Dyy=aT*qabs+(aL-aT)*qy2/qabs+Dmol*eps
                 	Dxy=(aL-aT)*(qx(itri))*(qy(itri))/qabs
                end if

                ! Face1
                qxface=S*qxval
           	 	qyface=qyval-ky(itri)*((rhorel-1)/12)*(5*con(k1)+5*con(k2)+2*con(k3))
            	dx= (x(k1)+x(k2)+x(k3))/3-(x(k1)+x(k2))/2
            	dy= (y(k1)+y(k2)+y(k3))/3-(y(k1)+y(k2))/2
            
            	face1_k1=(S**2*Dxx*Nx1+S*Dxy*Ny1)*dy-(Dyy*Ny1+S*Dxy*Nx1)*dx
            	face1_k2=(S**2*Dxx*Nx2+S*Dxy*Ny2)*dy-(Dyy*Ny2+S*Dxy*Nx2)*dx
            	face1_k3=(S**2*Dxx*Nx3+S*Dxy*Ny3)*dy-(Dyy*Ny3+S*Dxy*Nx3)*dx

            	qout=qxface*dy-qyface*dx;

            	if (qout >= 0) then
            		face1_k1=face1_k1-qout
            	else
            		face1_k2=face1_k2-qout
            	end if

            	! Face2
            	qxface=S*qxval
            	qyface=qyval-ky(itri)*((rhorel-1)/12)*(5*con(k1)+2*con(k2)+5*con(k3))
            
            	dx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2
            	dy= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2
            
            	face2_k1=(S**2*Dxx*Nx1+S*Dxy*Ny1)*dy-(Dyy*Ny1+S*Dxy*Nx1)*dx
            	face2_k2=(S**2*Dxx*Nx2+S*Dxy*Ny2)*dy-(Dyy*Ny2+S*Dxy*Nx2)*dx
            	face2_k3=(S**2*Dxx*Nx3+S*Dxy*Ny3)*dy-(Dyy*Ny3+S*Dxy*Nx3)*dx

            	qout=qxface*dely-qyface*delx

            	if (qout >= 0) then
            		face2_k1=face2_k1-qout
            	else
            		face2_k3=face2_k3-qout
            	end if

            	ap(k1)              = ap(k1) - face1_k1 - face2_k1
            	asup(k1,isup(k1))   = face1_k2 + face2_k2
            	asup(k1,isup(k1)+1) = face1_k3 + face2_k3
            	isup(k1)            = isup(k1) + 2
            end do
        end do



        conver=1;
    	conpre=connew;
    	do while (conver>tols)
    		do i = 1,N
    			do j = 1,maxsup
    				temp_sup(i,j) = connew(sup(i,j))
    			end do
    		end do
        	RHS=sum(asup*temp_sup,2);
        	connew=(eps*Volp*con+delt*RHS+BBs)/(eps*Volp+BCs+ap*delt);   
        	conver=maxval(abs(connew-conpre));
        	conpre=connew;
    	end do
    	con=connew

    end do


    print *, phi(1)
    print *, con(1)








	contains

	function interpolate1(x,y,z) result(res)
		implicit none
		real, dimension(1:2) :: x
		real, dimension(1:2) :: y
		integer, dimension(1:101) :: z
		real, dimension(1:101) :: res
		real :: slope

		integer  :: i

		
		slope = (y(2) - y(1))/(x(2) - x(1))
		

		do i = 1, size(z)
			res(i) = (z(i)-x(1))*slope + y(1)
		end do
	end function interpolate1





	






end program