subroutine FVM(L, x_0, N, gamma, Cv, Time, CFL, rho_l, u_l, p_l, rho_r, u_r, p_r)
	implicit none

	character(30), parameter					:: filename_output = "lastTimeStep.plt"

	integer										:: nt, N
	integer 									:: i, j, k
	integer, parameter							:: n_mon_point = 3, io_time_scan = 1488
	integer, dimension(n_mon_point)				:: i_mon_point, io_mon_point

	double precision							:: dt, t, Time, CFL, L, x_0

	double precision							:: gamma, Cv
	double precision							:: rho_l, u_l, p_l, H_l, c_l
	double precision							:: rho_r, u_r, p_r, H_r, c_r
	double precision 							:: rho_ll, u_ll, p_ll, rho_lr, u_lr, p_lr 			! for MUSCL
	double precision							:: rho_rl, u_rl, p_rl, rho_rr, u_rr, p_rr 			! for MUSCL
	double precision							:: dx, small_h, BIG_H
	double precision, dimension(N)			 	:: volume, volume_new
	double precision							:: x_left_piston, x_right_piston, amplitude_left, amplitude_right, &
												   frequency_left, frequency_right, u_left_piston, u_right_piston
	double precision							:: enthalpy, soundSpeed, veloPiston, xPiston
	double precision, dimension(3)				:: w, w_new, Flux_plus, Flux_minus
	double precision, dimension(0:N+1)			:: x, rho, u, p, H, rho_new, u_new, p_new			! x - centers of cells
	double precision, dimension(N+1)			:: surface_x, surface_height, u_dot
	double precision, dimension(n_mon_point)	:: x_mon_point
	
	double precision							:: maxu

	u_dot = 0.0
	small_h = 1.0
	BIG_H = 3.0 * small_h

	! creating grid
	dx = L / N
	call createQuasi1DGrid(N, L, x, dx, small_h, BIG_H, surface_x, surface_height, volume, volume_new)
	call writeChannelForm(N, x, surface_x, surface_height)
	! end creating grid

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
	! call BC_wall(rho(1), u(1), p(1), rho(0), u(0), p(0))
	call BC_piston(rho(1), u(1), p(1), rho(0), u(0), p(0), 0.0)

	! call BC_transmissive(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1))
	! call BC_wall(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1))
	call BC_piston(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1), 0.0)
	! end set initial conditions

	! set monitoring points
	x_mon_point(1) = 0.5
	x_mon_point(2) = 0.1
	x_mon_point(3) = 0.9
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

	! open file for time scan output
	open (io_time_scan, file = "time_scan.plt")
	write(io_time_scan,*) "variables = t, x, rho, u, p"
	! end open file for time scan output

	! set piston parameters
	x_left_piston = 0.0
	x_right_piston = L

	amplitude_left = 0.005
	amplitude_right = 0.005

	frequency_left = 500.0
	frequency_right = 500.0

	u_left_piston = 0.0
	u_right_piston = 0.0
	! end set piston parameters

	t = 0
	c_l = soundSpeed(gamma, rho(0), p(0), u(0))
	c_r = soundSpeed(gamma, rho(1), p(1), u(1))
	dt = CFL * dx / max(c_l, c_r)
	nt = Time / dt
	write(io_time_scan,*) "ZONE I=", nt, ", J=", N
	nt = 0
	do while (t < Time)

		t = t + dt
		nt = nt + 1

		u_left_piston = - veloPiston(amplitude_left, frequency_left, t) ! - veloPiston(amplitude_left, frequency_left, t)
		x_left_piston = - xPiston(amplitude_left, frequency_left, t) ! - xPiston(amplitude_left, frequency_left, t)
		u_right_piston = veloPiston(amplitude_right, frequency_right, t)
		x_right_piston = 1.0 + xPiston(amplitude_right, frequency_right, t)
		call deformationMesh(N, x, L, dx, dt, x_left_piston, x_right_piston, u_left_piston, u_right_piston, &
							 u_dot, surface_x, surface_height, small_h, BIG_H)

		do i = 1, N

			! first order
			! call roeScheme(gamma, Cv, rho(i-1), u(i-1), p(i-1), rho(i), u(i), p(i), u_dot(i), Flux_minus)
			! call roeScheme(gamma, Cv, rho(i), u(i), p(i), rho(i+1), u(i+1), p(i+1), u_dot(i+1), Flux_plus)
			call hllScheme(gamma, Cv, rho(i-1), u(i-1), p(i-1), rho(i), u(i), p(i), u_dot(i), Flux_minus)
			call hllScheme(gamma, Cv, rho(i), u(i), p(i), rho(i+1), u(i+1), p(i+1), u_dot(i+1), Flux_plus)
			! end first order 

			! second order
			! call TVD_MUSCL(i, N, rho, u, p, &
					 ! rho_ll, u_ll, p_ll, rho_lr, u_lr, p_lr, &
					 ! rho_rl, u_rl, p_rl, rho_rr, u_rr, p_rr)
			! ! write(*,*) rho_ll, u_ll, p_ll, rho_lr, u_lr, p_lr, rho_rl, u_rl, p_rl, rho_rr, u_rr, p_rr
			! ! call roeScheme(gamma, Cv, rho_ll, u_ll, p_ll, rho_lr, u_lr, p_lr, u_dot(i), Flux_minus)
			! ! call roeScheme(gamma, Cv, rho_rl, u_rl, p_rl, rho_rr, u_rr, p_rr, u_dot(i+1), Flux_plus)
			! call hllScheme(gamma, Cv, rho_ll, u_ll, p_ll, rho_lr, u_lr, p_lr, u_dot(i), Flux_minus)
			! call hllScheme(gamma, Cv, rho_rl, u_rl, p_rl, rho_rr, u_rr, p_rr, u_dot(i+1), Flux_plus)
			! end second order

			H(i) = enthalpy(gamma, rho(i), p(i), u(i))
			call calcVariables(rho(i), u(i), p(i), H(i), w)

			volume_new(i) = (surface_height(i) + surface_height(i+1)) / 2.0 * (surface_x(i+1) - surface_x(i))
			! write(*,*) x(i), volume(i), volume_new(i)

			call eulerTimeScheme(w_new, w, volume(i), volume_new(i), dt, Flux_plus, Flux_minus, &
								 surface_height(i), surface_height(i+1), p(i-1), p(i+1))

			call calcVariablesFromVector(gamma, rho_new(i), u_new(i), p_new(i), w_new)

		end do

		rho = rho_new
		u = u_new
		p = p_new
		volume = volume_new

		call findMamaCell(N, x, surface_x, n_mon_point, x_mon_point, i_mon_point)
		do k = 1, n_mon_point
			write(io_mon_point(k),*) t, rho(i_mon_point(k)), u(i_mon_point(k)), p(i_mon_point(k))
		end do

		! set boundary conditions
		! call BC_transmissive(rho(1), u(1), p(1), rho(0), u(0), p(0))
		! call BC_wall(rho(1), u(1), p(1), rho(0), u(0), p(0))
		call BC_piston(rho(1), u(1), p(1), rho(0), u(0), p(0), u_left_piston)

		! call BC_transmissive(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1))
		! call BC_wall(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1))
		call BC_piston(rho(N), u(N), p(N), rho(N+1), u(N+1), p(N+1), u_right_piston)
		! end set boundary conditions

		if (mod(nt, 1000) == 0) then
			write(*,*) nt, t
			write(*,*) x_left_piston, x_right_piston - 1.0
			write(*,*) u_left_piston, u_right_piston
			write(*,*) rho(N/2), u(N/2), p(N/2)
		end if
		
		! call writeTimeScan(t, N, x, rho, u, p)

	end do
	close(io_time_scan)

	write(*,*) "Time = ", t
	write(*,*) u_left_piston, x_left_piston
	open(685, file = "udots.plt")
	write(685,*) "variables = x, u"
	do i = 1, N + 1
		write(685,*) surface_x(i), u_dot(i)
	end do
	call writeField(filename_output, N, x, surface_x, rho, u, p)

end subroutine

! Euler explicit scheme
subroutine eulerTimeScheme(w_new, w, volume, volume_new, dt, Flux_plus, Flux_minus, surface_height_l, surface_height_r, p_l, p_r)
	implicit none

	double precision					:: p_l, p_r, dt, volume, volume_new
	double precision					:: surface_height_l, surface_height_r
	double precision, dimension(3)		:: w, w_new, Flux_minus, Flux_plus

	w_new = w * volume / volume_new - dt / volume_new * (Flux_plus * surface_height_r - Flux_minus * surface_height_l)
	w_new(2) = w_new(2) + dt / volume_new * (p_l + p_r) / 2.0 * (surface_height_r - surface_height_l)

end subroutine