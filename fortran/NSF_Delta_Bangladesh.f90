program main

use utils

implicit none

real*4 :: x1 = 1
real*8 :: x2 = 10
real*8 :: array(10)

real*8 :: x(5,10)

integer :: i,j

call linspace(x1,x2,array)

do i = 1,10
	print *, array(i)
end do

call read_in()



end program main


! subroutine read_in
! 	implicit none

! 	real*8 :: Breath = 300000
! 	integer :: cols = 101

! 	real*8 :: Delx, Dely
! 	real*8, allocatable, dimension(:) :: xcols

! 	real*8 :: AR = 80

! 	real*8 :: delt = 100 ! timstep days

! 	real*8 :: consea =1 ! salt concentration in sea
! 	real*8 :: rhorel = 1.025 ! relative density of saturated saline

! 	real*8 :: eps = 0.35 !aquifer porosity
! 	real*8 :: sto = 0.001 !aquifer storage

! 	! hydraulic conductivity m/day
! 	real*8 :: kxval = 10
! 	real*8 :: kyval 
! 	! diffusion and dispersion coefficients
! 	real*8 :: Dmol = 0.00001 !m2/day
! 	real*8 :: aL = 50
! 	real*8 :: aT 

! 	!USER DEFINED SYNTHETIC STRATIGRAPHY
! 	real*8 :: sealevel = 300 ! sea level setting
! 	real*8 :: basement = 250 ! depth of the basement

! 	real*8 :: foreslope = 0.05 !foreset slope
! 	real*8 :: topslope = 0.0005 ! fluvial topset slope

! 	real*8 :: xtoein 
! 	real*8 :: xtoe 

! 	real*8 :: toevel = 0.04 ! toe velocity over basement m/day
! 	real*8 :: xshore

! 	real*8, allocatable, dimension(:,:) :: eta
! 	real*8, allocatable, dimension(:,:) :: etabot

! 	real*8, dimension(1,101) :: aa
	

! 	integer :: tstep = 1
! 	integer :: totstep

! 	Delx = Breath/(cols-1)
! 	Dely = Delx/AR

! 	kyval = kxval/100
! 	aT = aL/10

! 	xtoein = (sealevel - basement)/foreslope
! 	xtoe = xtoein

! 	! do while (xtoe < 0.8*Breath)
! 	! 	xtoe = xtoein + toevel * (tstep-1)*delt
! 	! 	xshore = xtoe - (sealevel - basement)/foreslope
! 	! 	eta(tstep, :) = min ( max((xshore-xcols)*topslope+sealevel,sealevel),&
!  !                       max((xtoe-xcols)*foreslope+basement,basement) )
! 	! 	etabot(tstep,:)= zeros(1,cols)
!  !    	tstep=tstep+1
!  !    end do
!  !    totstep = tstep - 1
!  	aa = zeros(1,101)
! end subroutine





	






