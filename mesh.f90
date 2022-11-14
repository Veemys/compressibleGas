! create 1D grid
subroutine create1DGrid(N, x, dx)
	implicit none

	integer									:: N
	integer 								:: i
	double precision						:: dx
	double precision, dimension(0:N + 1)	:: x

	x(0) = - dx / 2.0
	do i = 1, N + 1

		x(i) = x(i-1) + dx

	end do

end subroutine

! create quasi 1D grid
!
!	--------------------
!	|                  |
!	|				   |
!	i         i       i+1
!	|                  |
!	|                  |
!	--------------------
!
subroutine createQuasi1DGrid(N, x, dx, small_h, BIG_H, surface_x)
	implicit none

	integer									:: N
	integer 								:: i
	double precision						:: dx, small_h, BIG_H
	double precision, dimension(0:N + 1)	:: x
	double precision, dimension(N + 1)		:: surface_x, height

	x(0) = - dx / 2.0
	surface_x(1) = 0.0
	do i = 1, N + 1

		! calculation of the cells centers
		x(i) = x(i-1) + dx
		
		! calculation of the surfaces position
		surface_x(i) = x(i-1) + dx / 2.0

		! calculation of the height of channel
		if (surface_x(i) < 1.0) then
		
			height(i) = small_h
		
		else if (surface_x(i) > 2.0) then
		
			height(i) = BIG_H
		
		else
		
			height(i) = ! осталось придумать прикольную функцию
		
		end if

	end do

end subroutine

! deformation quasi 1D grid
subroutine deformationMesh(N, L, dt, u_left_piston, u_right_piston, u_surface)
	implicit none
	
	integer								:: i, N
	double precision					:: L, dt
	double precision					:: u_left_piston, u_right_piston
	double precision, dimension(N+1)	:: u_surface, surface_x			! u_surface - sounds of surfaces
	double precision, dimension(0:N+1)	::
	
	! calculation of sounds of surfaces and of the surfaces coordinates
	do i = 1, N + 1
		
		u_surface(i) = u_left_piston + (u_right_piston - u_left_piston) / L * surface_x(i)
		surface_x(i) = surface_x(i) + u_surface(i) * dt
		
	end do
	
	! calculation of the centers of cells
	x(0) = - dx / 2.0
	do i = 1, N
	
		x(i) = (surface_x(i) + surface_x(i+1)) / 2.0
		
	end do
	L = surface_x(N+1)
	x(N) = L + dx / 2.0
	
end subroutine