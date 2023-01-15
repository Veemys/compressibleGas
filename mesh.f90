! create 1D grid
subroutine create1DGrid(N, x, dx)
	implicit none

	integer									:: N
	integer 								:: i
	double precision						:: dx
	double precision, dimension(0:N+1)	:: x

	x(0) = - dx / 2.0
	do i = 1, N + 1

		x(i) = x(i-1) + dx

	end do

end subroutine

! create quasi 1D grid
!
!	---------------------
!	|                   |
!	|				    |
!	i         i        i+1				<- cell and her surfaces
!	|                   |
!	|                   |
!	---------------------
!
subroutine createQuasi1DGrid(N, L, x, dx, small_h, BIG_H, surface_x, surface_height, volume, volume_new)
	implicit none

	integer									:: N
	integer 								:: i
	double precision						:: dx, small_h, BIG_H, L
	double precision, dimension(0:N+1)		:: x
	double precision, dimension(N+1)		:: surface_x, surface_height
	double precision, dimension(N)			:: volume, volume_new

	x(0) = - dx / 2.0
	surface_x(1) = 0.0
	do i = 1, N + 1

		! calculation of the cells centers
		x(i) = x(i-1) + dx
		
		! calculation of the surfaces position
		surface_x(i) = x(i-1) + dx / 2.0

		! calculation of the height of channel
		if (surface_x(i) < L / 3.0) then
		
			surface_height(i) = small_h
		
		else if (surface_x(i) > 2.0 * L / 3.0) then
		
			surface_height(i) = BIG_H
		
		else
		
			surface_height(i) = small_h + (BIG_H - small_h) * sin(surface_x(i) - (L / 3.0)) / sin(L / 3.0)
		
		end if

	end do
	
	volume = 0.0
	volume_new = 0.0
	do i = 1, N
		
		volume(i) = (surface_height(i) + surface_height(i+1)) / 2.0 * (surface_x(i+1) - surface_x(i))
		volume_new(i) = (surface_height(i) + surface_height(i+1)) / 2.0 * (surface_x(i+1) - surface_x(i))
		
	end do

end subroutine

! deformation quasi 1D grid
subroutine deformationMesh(N, x, L, dx, dt, x_left_piston, x_right_piston, u_left_piston, u_right_piston, &
						   u_surface, surface_x, surface_height, small_h, BIG_H)
	implicit none

	integer								:: i, N
	double precision					:: L, dt, dx
	double precision					:: x_left_piston, x_right_piston, u_left_piston, u_right_piston
	double precision					:: small_h, BIG_H
	double precision, dimension(N+1)	:: u_surface, surface_x, surface_height			! u_surface - sounds of surfaces
	double precision, dimension(0:N+1)	:: x

	! calculation of sounds of surfaces and of the surfaces coordinates
	surface_x(1) = x_left_piston
	u_surface(1) = u_left_piston
	surface_x(N+1) = x_right_piston
	u_surface(N+1) = u_right_piston
	do i = 2, N

		! u_surface(i) = u_surface(1) + (u_right_piston - u_left_piston) / L * surface_x(i)
		! u_surface(i) = u_surface(1) + ((surface_x(i) - surface_x(1)) * (u_surface(N+1) - u_surface(1))) / (surface_x(N+1) - surface_x(1))
		u_surface(i) = u_surface(1) + (u_surface(N+1) - u_surface(1)) / (surface_x(N+1) - surface_x(1)) * surface_x(i)
		surface_x(i) = surface_x(i) + u_surface(i) * dt

	end do

	! calculation of the height of channel
	do i = 1, N + 1

		if (surface_x(i) < 1.0 / 3.0) then

			surface_height(i) = small_h

		else if (surface_x(i) > 2.0 / 3.0) then

			surface_height(i) = BIG_H

		else

			surface_height(i) = small_h + (BIG_H - small_h) * sin(surface_x(i) - (L / 3.0)) / sin(L / 3.0)

		end if

	end do

	! calculation of the centers of cells
	x(0) = surface_x(1) - dx / 2.0
	do i = 1, N

		x(i) = (surface_x(i) + surface_x(i+1)) / 2.0

	end do
	L = surface_x(N+1) - surface_x(1)
	x(N+1) = surface_x(N+1) + dx / 2.0 ! L + dx / 2.0

end subroutine