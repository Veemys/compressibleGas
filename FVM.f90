subroutine FVM(L, x_0, N, gamma, Cv, Time, CFL, rho_l, u_l, p_l, rho_r, u_r, p_r)
	implicit none

	character(30), parameter					:: filename_output = "lastTimeStep.plt"

	integer										:: nt, N
	integer 									:: i, j, k
	integer, parameter							:: n_mon_point = 2
	integer, dimension(n_mon_point)				:: i_mon_point, io_mon_point

	double precision							:: dt, t, Time, CFL, L, x_0

	double precision							:: gamma, Cv
	double precision							:: rho_l, u_l, p_l, H_l, c_l
	double precision							:: rho_r, u_r, p_r, H_r, c_r
	double precision							:: dx, volume, volume_new, small_h, BIG_H
	double precision							:: x_left_piston, x_right_piston, amplitude_left, amplitude_right, &
												   frequency_left, frequency_right, u_left_piston, u_right_piston
	double precision							:: enthalpy, soundSpeed, veloPiston, xPiston
	double precision, dimension(3)				:: w, w_new, Flux_plus, Flux_minus
	double precision, dimension(0:N + 1)		:: x, rho, u, p, H, rho_new, u_new, p_new			! x - centers of cells
	double precision, dimension(N + 1)			:: surface_x, surface_height, u_dot
	double precision, dimension(n_mon_point)	:: x_mon_point

	u_dot = 0.0
	small_h = 1.0
	BIG_H = small_h

	! creating grid
	dx = L / N
	call createQuasi1DGrid(N, x, dx, small_h, BIG_H, surface_x, surface_height)
	call writeChannelForm(N, surface_x, surface_height)

	! set initial conditions
	do i = 1, N

		if (x(i) <= x_0) then
			rho(i) = rho_l
			u(i) = u_l
			p(i) = p_l
		else if (x(i) > x_0) then
			rho(i) = rho_r
			u(i) = u_r
			p(i) = p_r
		end if

	end do

	! call BC_transmissive(rho(1), u(1), p(1), rho(0), u(0), p(0))
	call BC_wall(rho(1), u(1), p(1), rho(0), u(0), p(0))
	! call BC_piston(rho(1), u(1), p(1), rho(0), u(0), p(0), u_p)
	
	! call BC_transmissive(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1))
	call BC_wall(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1))
	! call BC_piston(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1), u_p)
	! end set initial conditions
	
	! set monitoring points
	x_mon_point(1) = 0.3
	x_mon_point(2) = 0.6
	io_mon_point = 6969
	call findMamaCell(N, x, surface_x, n_mon_point, x_mon_point, i_mon_point)
	do i = 1, n_mon_point
		
		io_mon_point(i) = io_mon_point(i) + i
		if (i == 1) then
			open(io_mon_point(i), file = "monitoring_points_1.plt")
		else if (i == 2) then
			open(io_mon_point(i), file = "monitoring_points_2.plt")
		else if (i == 3) then
			open(io_mon_point(i), file = "monitoring_points_3.plt")
		else if (i == 4) then
			open(io_mon_point(i), file = "monitoring_points_4.plt")
		else if (i == 5) then
			open(io_mon_point(i), file = "monitoring_points_5.plt")
		end if
		write(io_mon_point(i),*) "variables = t, rho, u, p"
		
	end do
	! end set monitoring points

	! set piston parameters
	x_left_piston = 0.0
	x_right_piston = L
	
	amplitude_left = 0.1
	amplitude_right = 0.0
	
	frequency_left = 1073
	frequency_right = 0.0
	
	u_left_piston = 0.0
	u_right_piston = 0.0
	! end set piston parameters

	t = 0
	nt = 0
	c_l = soundSpeed(gamma, rho(0), p(0), u(0))
	c_r = soundSpeed(gamma, rho(1), p(1), u(1))
	dt = CFL * dx / max(c_l, c_r)
	do while (t < Time)
		
		t = t + dt
		nt = nt + 1
		
		u_left_piston = veloPiston(amplitude_left, frequency_left, t)
		x_left_piston = xPiston(x_left_piston, u_left_piston, dt)
		call deformationMesh(N, x, L, dx, dt, u_left_piston, u_right_piston, u_dot, surface_x)
		
		do i = 1, N
			
			! call roeScheme(gamma, Cv, rho(i-1), u(i-1), p(i-1), rho(i), u(i), p(i), u_dot(i), Flux_minus)
			! call roeScheme(gamma, Cv, rho(i), u(i), p(i), rho(i+1), u(i+1), p(i+1), u_dot(i+1), Flux_plus)
			call hllScheme(gamma, Cv, rho(i-1), u(i-1), p(i-1), rho(i), u(i), p(i), u_dot(i), Flux_minus)
			call hllScheme(gamma, Cv, rho(i), u(i), p(i), rho(i+1), u(i+1), p(i+1), u_dot(i+1), Flux_plus)

			! c_l = soundSpeed(gamma, (rho(i-1) + rho(i)) / 2.0, (p(i-1) + p(i)) / 2.0, (u(i-1) + u(i)) / 2.0)
			! c_r = soundSpeed(gamma, (rho(i+1) + rho(i)) / 2.0, (p(i+1) + p(i)) / 2.0, (u(i+1) + u(i)) / 2.0)
			! dt = CFL * (surface_x(i+1) - surface_x(i)) / max(c_l, c_r)

			H(i) = enthalpy(gamma, rho(i), p(i), u(i))
			call calcVariables(rho(i), u(i), p(i), H(i), w)
			
			volume = (surface_height(i) + surface_height(i+1)) / 2.0 * (surface_x(i+1) - surface_x(i))
			volume_new = (surface_height(i) + surface_height(i+1)) / 2.0 * (surface_x(i+1) - surface_x(i))

			call eulerTimeScheme(w_new, w, volume, volume_new, dt, Flux_plus, Flux_minus, surface_height(i), &
								 surface_height(i + 1), p(i-1), p(i+1))

			call calcVariablesFromVector(gamma, rho_new(i), u_new(i), p_new(i), w_new)

		end do

		rho = rho_new
		u = u_new
		p = p_new
		
		call findMamaCell(N, x, surface_x, n_mon_point, x_mon_point, i_mon_point)
		do k = 1, n_mon_point
			write(io_mon_point(k),*) t, rho(i_mon_point(k)), u(i_mon_point(k)), p(i_mon_point(k))
		end do
		
		! set boundary conditions
		
		! call BC_transmissive(rho(1), u(1), p(1), rho(0), u(0), p(0))
		call BC_wall(rho(1), u(1), p(1), rho(0), u(0), p(0))
		! call BC_piston(rho(1), u(1), p(1), rho(0), u(0), p(0), u_p)
		
		! call BC_transmissive(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1))
		call BC_wall(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1))
		! call BC_piston(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1), u_p)
		
		! end set boundary conditions
		
		write(*,*) nt, t

	end do
	
	write(*,*) "Time = ", t

	call writeField(filename_output, N, x, rho, u, p)

end subroutine

! Euler explicit scheme
subroutine eulerTimeScheme(w_new, w, volume, volume_new, dt, Flux_plus, Flux_minus, surface_height_l, surface_height_r, p_l, p_r)
	implicit none
	
	double precision					:: p_l, p_r, dt, volume, volume_new
	double precision					:: surface_height_l, surface_height_r
	double precision, dimension(3)		:: w, w_new, Flux_minus, Flux_plus
	
	w_new = w * volume / volume_new - dt / volume_new * (Flux_plus * surface_height_r - Flux_minus * surface_height_l)
	w_new(2) = w_new(2) + (p_l + p_r) / 2.0 * (surface_height_r - surface_height_l)

end subroutine

! implicit predictor-corrector time scheme
subroutine implicitTimeScheme(N, gamma, rho, u, p, w_new, w, volume, volume_new, dt, Flux_plus, Flux_minus, &
							  surface_height_l, surface_height_r, p_l, p_r)
	implicit none
	
	integer								:: N
	double precision					:: gamma, p_l, p_r, dt, volume, volume_new
	double precision					:: surface_height_l, surface_height_r
	double precision, dimension(3)		:: w, w_new, Flux_minus, Flux_plus
	double precision, dimension(0:N+1)	:: rho, u, p
	
	call calcVariablesFromVector(gamma, rho, u, p, w)
	call predictor(N)
	call corrector(N)

end subroutine

! predictor
subroutine predictor(N)
	implicit none
	
	integer								:: N
	double precision					:: gamma, p_l, p_r, dt, volume, volume_new
	double precision					:: surface_height_l, surface_height_r
	double precision, dimension(3)		:: w, w_new, Flux_minus, Flux_plus
	double precision, dimension(0:N+1)	:: rho, u, p
	
	

end subroutine

! corrector
subroutine corrector(N)
	implicit none
	
	integer								:: N
	double precision					:: gamma, p_l, p_r, dt, volume, volume_new
	double precision					:: surface_height_l, surface_height_r
	double precision, dimension(3)		:: w, w_new, Flux_minus, Flux_plus
	double precision, dimension(0:N+1)	:: rho, u, p
	
	

end subroutine