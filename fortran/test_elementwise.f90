program test
	implicit none

	real, dimension(3) :: x = (/1,2,3/), y = (/4,5,6/)
	real, dimension(3) :: z
	real, dimension(3,3) :: t
	z = x*y
	t = reshape((/ 1, 2, 3, 2, 3, 1, 3, 1, 2 /), shape(t))
	! print *, z
	! print *, sum(z)
	print *, t
	print *, sum(t,2)

end program test