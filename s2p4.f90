program s2p4

	implicit none

	integer, parameter :: nx = 1000, ny = 100, nz = 4
	integer :: i,j,k
	double precision, dimension(0:nx,0:ny,0:nz) :: u, T, x, y, z
		
	open(unit=20, file="profile.dat")
        read(20,*)
        read(20,*)

	do k = 0, nz
		do j = 0, ny
			do i = 0, nx
				read(20, '(5E30.18)') x(i,j,k), y(i,j,k), z(i,j,k), u(i,j,k), T(i,j,k)
			end do
		end do
	end do

	open(unit=15, file="Turb_inflow_bl00.dat")
	open(unit=13, file="Turb_inflow_bl01.dat")
	open(unit=11, file="Turb_inflow_bl02.dat")
	open(unit=10, file="Turb_inflow_bl03.dat")

	do k = 0, nz
		do j = 0, ny
			do i = 0, nx
				if ( j .le. ny/2 .and. k .le. nz/2 ) then
					write(15,'(5E30.18)') x(i,j,k), y(i,j,k), z(i,j,k), u(i,j,k), T(i,j,k)
				elseif ( j .le. ny/2 .and. k .gt. nz/2 ) then
					write(13, '(5E30.18)') x(i,j,k), y(i,j,k), z(i,j,k), u(i,j,k), T(i,j,k)
				elseif ( j .gt. ny/2  .and. k .le. nz/2 ) then
					write(11, '(5E30.18)') x(i,j,k), y(i,j,k), z(i,j,k), u(i,j,k), T(i,j,k)
				elseif ( j .gt. ny/2 .and. k .gt. nz/2 ) then
					write(10, '(5E30.18)') x(i,j,k), y(i,j,k), z(i,j,k), u(i,j,k), T(i,j,k)
				end if
			end do
		end do
	end do

	close(20)	
	close(15)
	close(13)
	close(11)
	close(10)

end program	

