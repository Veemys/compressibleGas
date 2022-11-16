subroutine FVM(L, x_0, N, gamma, Cv, Time, CFL, rho_l, u_l, p_l, rho_r, u_r, p_r)
	implicit none

	character(30), parameter				:: filename_output = "lastTimeStep.plt"

	integer									:: nt, N
	integer 								:: i, j, k

	double precision						:: dt, t, Time, CFL, L, x_0

	double precision						:: gamma, Cv
	double precision						:: rho_l, u_l, p_l, H_l, c_l
	double precision						:: rho_r, u_r, p_r, H_r, c_r
	double precision						:: dx, volume, volume_new, small_h, BIG_H
	double precision						:: enthalpy, soundSpeed
	double precision, dimension(3)			:: w, w_new, Flux_plus, Flux_minus
	double precision, dimension(0:N + 1)	:: x, rho, u, p, H, rho_new, u_new, p_new			! x - centers of cells
	double precision, dimension(N + 1)		:: surface_x, surface_height, u_dot

	u_dot = 0.0
	small_h = 1.0
	BIG_H = small_h

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

	rho(0) = rho(1)
	u(0) = u(1)
	p(0) = p(1)

	rho(N+1) = rho(N)
	u(N+1) = u(N)
	p(N+1) = p(N)
	! end set initial conditions

	t = 0
	nt = 0
	c_l = soundSpeed(gamma, rho(0), p(0), u(0))
	c_r = soundSpeed(gamma, rho(1), p(1), u(1))
	dt = CFL * dx / max(c_l, c_r)
	do while (t < Time)
		
		t = t + dt
		nt = nt + 1
		
		do i = 1, N
			
			! call roeScheme(gamma, Cv, rho(i-1), u(i-1), p(i-1), rho(i), u(i), p(i), u_dot(i), Flux_minus)
			! call roeScheme(gamma, Cv, rho(i), u(i), p(i), rho(i+1), u(i+1), p(i+1), u_dot(i+1), Flux_plus)
			call hllScheme(gamma, Cv, rho(i-1), u(i-1), p(i-1), rho(i), u(i), p(i), u_dot, Flux_minus)
			call hllScheme(gamma, Cv, rho(i), u(i), p(i), rho(i+1), u(i+1), p(i+1), u_dot, Flux_plus)

			c_l = soundSpeed(gamma, (rho(i-1) + rho(i)) / 2.0, (p(i-1) + p(i)) / 2.0, (u(i-1) + u(i)) / 2.0)
			c_r = soundSpeed(gamma, (rho(i+1) + rho(i)) / 2.0, (p(i+1) + p(i)) / 2.0, (u(i+1) + u(i)) / 2.0)
			! dt = CFL * (surface_x(i+1) - surface_x(i)) / max(c_l, c_r)

			H(i) = enthalpy(gamma, rho(i), p(i), u(i))
			call calcVariables(rho(i), u(i), p(i), H(i), w)
			
			volume = (surface_height(i) + surface_height(i+1)) / 2.0 * (surface_x(i+1) - surface_x(i))
			volume_new = (surface_height(i) + surface_height(i+1)) / 2.0 * (surface_x(i+1) - surface_x(i))

			call eulerTimeScheme(w_new, w, volume, volume_new, dt, Flux_plus, Flux_minus, surface_height(i), &
								 surface_height(i + 1), p(i-1), p(i+1))
			! call implicitTimeScheme(w_new, w, volume, volume_new, dt, Flux_plus, Flux_minus, surface_height(i), &
			! 						  surface_height(i + 1), p(i-1), p(i+1))

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

! implicit predictor-corrector scheme
subroutine implicitTimeScheme(N, gamma, rho, u, p, w_new, w, volume, volume_new, dt, Flux_plus, Flux_minus, surface_height_l, surface_height_r, p_l, p_r)
	implicit none
	
	integer								:: N
	double precision					:: gamma, p_l, p_r, dt, volume, volume_new
	double precision					:: surface_height_l, surface_height_r
	double precision, dimension(3)		:: w, w_new, Flux_minus, Flux_plus
	double precision, dimension(0:N+1)	:: rho, u, p
	
	call calcVariablesFromVector(gamma, rho, u, p, w)
	call predictor()
	call corrector()

end subroutine

! predictor
subroutine predictor()
	implicit none
	
	integer								:: N
	double precision					:: gamma, p_l, p_r, dt, volume, volume_new
	double precision					:: surface_height_l, surface_height_r
	double precision, dimension(3)		:: w, w_new, Flux_minus, Flux_plus
	double precision, dimension(0:N+1)	:: rho, u, p
	
	

end subroutine

! corrector
subroutine corrector()
	implicit none
	
	integer								:: N
	double precision					:: gamma, p_l, p_r, dt, volume, volume_new
	double precision					:: surface_height_l, surface_height_r
	double precision, dimension(3)		:: w, w_new, Flux_minus, Flux_plus
	double precision, dimension(0:N+1)	:: rho, u, p
	
	

end subroutine