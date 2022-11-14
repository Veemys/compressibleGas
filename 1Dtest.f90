subroutine test1D(L, x_0, N, gamma, Cv, Time, CFL, rho_l, u_l, p_l, rho_r, u_r, p_r)
	implicit none

	character(30), parameter				:: filename = "lastTimeStep.plt"

	integer									:: nt, N
	integer 								:: i, j, k

	double precision						:: dx, dt, t, Time, CFL, L, x_0

	double precision						:: gamma, Cv
	double precision						:: rho_l, u_l, p_l, H_l, c_l
	double precision						:: rho_r, u_r, p_r, H_r, c_r
	double precision						:: u_dot
	double precision						:: enthalpy, soundSpeed
	double precision, dimension(3)			:: w, w_new, Flux_plus, Flux_minus
	double precision, dimension(0:N + 1)	:: x, rho, u, p, H, rho_new, u_new, p_new

	u_dot = 0.0

	dx = L / N
	call create1DGrid(N, x, dx)

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

	rho(0) = rho(1)
	u(0) = u(1)
	p(0) = p(1)

	rho(N+1) = rho(N)
	u(N+1) = u(N)
	p(N+1) = p(N)
	! end set initial conditions

	nt = 0
	t = 0
	c_l = soundSpeed(gamma, rho(0), p(0), u(0))
	c_r = soundSpeed(gamma, rho(1), p(1), u(1))
	dt = CFL * dx / max(c_l, c_r) / 100.0
	do while (t < Time)
		
		t = t + dt
		nt = nt + 1
		
		do i = 1, N
			
			call roeScheme(gamma, Cv, rho(i-1), u(i-1), p(i-1), rho(i), u(i), p(i), u_dot, Flux_minus)
			call roeScheme(gamma, Cv, rho(i), u(i), p(i), rho(i+1), u(i+1), p(i+1), u_dot, Flux_plus)
			! call hllScheme(gamma, Cv, rho(i-1), u(i-1), p(i-1), rho(i), u(i), p(i), u_dot, Flux_minus)
			! call hllScheme(gamma, Cv, rho(i), u(i), p(i), rho(i+1), u(i+1), p(i+1), u_dot, Flux_plus)

			c_l = soundSpeed(gamma, (rho(i-1) + rho(i)) / 2.0, (p(i-1) + p(i)) / 2.0, (u(i-1) + u(i)) / 2.0)
			c_r = soundSpeed(gamma, (rho(i+1) + rho(i)) / 2.0, (p(i+1) + p(i)) / 2.0, (u(i+1) + u(i)) / 2.0)
			dt = CFL * dx / max(c_l, c_r) / 100.0

			H(i) = enthalpy(gamma, rho(i), p(i), u(i))
			call calcVariables(rho(i), u(i), p(i), H(i), w)

			w_new = w - dt / dx * (Flux_plus - Flux_minus)
			call calcVariablesFromVector(gamma, rho_new(i), u_new(i), p_new(i), w_new)

		end do

		rho = rho_new
		u = u_new
		p = p_new
		
		! set boundary conditions
		rho(0) = rho(1)
		u(0) = u(1)
		p(0) = p(1)

		rho(N+1) = rho(N)
		u(N+1) = u(N)
		p(N+1) = p(N)
		! end set boundary conditions
		
		write(*,*) nt, t

	end do

	write(*,*) "Time = ", t

	call writeField(filename, N, x, rho, u, p)

end subroutine