module utils

implicit none

real*8, public :: Breath = 300000
integer, public :: cols = 101
real*8, public :: Delx, Dely
real*8, public, allocatable, dimension(:) :: xcols
real*8, public :: AR = 80
real*8, public :: delt = 100 ! timstep days
real*8, public :: consea =1 ! salt concentration in sea
real*8, public :: rhorel = 1.025 ! relative density of saturated saline
real*8, public :: eps = 0.35 !aquifer porosity
real*8, public :: sto = 0.001 !aquifer storage
! hydraulic conductivity m/day
real*8, public :: kxval = 10
real*8, public :: kyval 
! diffusion and dispersion coefficients
real*8, public :: Dmol = 0.00001 !m2/day
real*8, public :: aL = 50
real*8, public :: aT 
!USER DEFINED SYNTHETIC STRATIGRAPHY
real*8, public :: sealevel = 300 ! sea level setting
real*8, public :: basement = 250 ! depth of the basement
real*8, public :: foreslope = 0.05 !foreset slope
real*8, public :: topslope = 0.0005 ! fluvial topset slope
real*8, public :: xtoein 
real*8, public :: xtoe 
real*8 , public :: toevel = 0.04 ! toe velocity over basement m/day
real*8 , public:: xshore
real*8, public, allocatable, dimension(:,:) :: eta
real*8, public, allocatable, dimension(:,:) :: temp
real*8, public, allocatable, dimension(:,:) :: etabot

integer, public :: tstep = 1
integer, public :: totstep


contains

function zeros(m,n) result(array)
	integer :: m, n, i, j
	real*8, dimension(m,n) :: array
	do i = 1,m
		do j = 1,n
			array(i,j) = 0
		end do
	end do
end function

function interpolate1(x,y,z) result(res)
	real*8, dimension(1:2) :: x
	real*8, dimension(1:2) :: y
	integer, dimension(1:101) :: z
	real*8, dimension(1:101) :: res
	real*8 :: slope

	integer  :: i
	
	slope = (y(2) - y(1))/(x(2) - x(1))
	
	do i = 1, size(z)
		res(i) = (z(i)-x(1))*slope + y(1)
	end do
end function interpolate1


subroutine linspace(x1, x2, array)
	! implicit none
	real*4, intent(in) :: x1
	real*8, intent(in) :: x2
	real*8, intent(out) :: array(:)

	real*8 :: range
	integer :: n, i
	n = size(array)

	range = x2 - x1

	if (n == 0) then
		return
	end if

	if (n == 1) then
		array(1) = x1
		return
	end if

	do i = 1,n
		array(i) = x1 + range*(i-1)/(n-1)
	end do
end subroutine

subroutine initialize(array, value)
	real*8, allocatable, dimension(:,:), intent(out) :: array
	integer :: m, n, i, j
	real*8 :: value
	m = size(array,1)
	n = size(array,2)
	do i = 1, m
		do j = 1, n
		 array(i,j) = value
		end do
	end do
end subroutine

subroutine read_in()

	Delx = Breath/(cols-1)
	Dely = Delx/AR
	kyval = kxval/100
	aT = aL/10
	xtoein = (sealevel - basement)/foreslope
	xtoe = xtoein

	allocate(xcols(cols))
	call linspace(0.0, breath, xcols)
	allocate(eta(tstep,cols))
	allocate(temp(1,cols))
	! print *, shape(eta)
	! print *, max((xshore-xcols)*topslope+sealevel,sealevel)
	do while (xtoe < 0.8*Breath)
		xtoe = xtoein + toevel * (tstep-1)*delt
		xshore = xtoe - (sealevel - basement)/foreslope
		! temp = min(max((xshore-xcols)*topslope+sealevel,sealevel),&
                       max((xtoe-xcols)*foreslope+basement,basement))
		! temp = reshape(min(max((xshore-xcols)*topslope+sealevel,sealevel),&
  !                      max((xtoe-xcols)*foreslope+basement,basement)), (/1,101/))
		! ! etabot(tstep,:)= zeros(1,cols)
		! eta(tstep,:) = temp
		
		if (allocated(eta)) then
        	deallocate(eta)
    	endif
    	allocate(eta(tstep, cols))
    	! call move_alloc(temp, eta(tstep,cols))
    	eta(tstep,:) = min(max((xshore-xcols)*topslope+sealevel,sealevel),&
                       max((xtoe-xcols)*foreslope+basement,basement))
    	tstep=tstep+1
    	
    end do
    totstep = tstep - 1

    print *, shape(eta)
    print *, eta(2,:)
end subroutine

end module




