subroutine test1D()
implicit none

character(30), parameter				:: filename = "lastTimeStep.plt"

integer,parameter						:: N = 100
integer 								:: i, j, k

double precision						:: dx, dt, t, Time, CFL, L

double precision						:: gamma, Cv
double precision						:: rho_l, u_l, H_l, p_l, c_l
double precision						:: rho_r, u_r, H_r, p_r, c_r
double precision						:: u_dot
double precision						:: enthalpy, soundSpeed
double precision, dimension(3)			:: w, w_new, Flux_plus, Flux_minus
double precision, dimension(0:N + 1)	:: x, rho, u, p, H, rho_new, u_new, p_new

Cv = 1005.
gamma = 1.4
u_dot = 0.0
L = 1
Time = 5e-4
CFL = 0.1

dx = L / N
call create1DGrid(N, x, dx)

! set initial conditions
do i = 1, N

	if (x(i) <= L / 2.0) then
		rho(i) = 1.2
		u(i) = 0.
		p(i) = 1e5
	else if (x(i) > L / 2.0) then
		rho(i) = 0.047
		u(i) = 0.
		p(i) = 4000.
	end if

end do

rho(0) = rho(1)
u(0) = - u(1)
p(0) = p(1)

rho(N+1) = rho(N)
u(N+1) = - u(N)
p(N+1) = p(N)
! end set initial conditions

t = 0
c_l = soundSpeed(gamma, rho(0), p(0), u(0))
c_r = soundSpeed(gamma, rho(1), p(1), u(1))
dt = CFL * dx / max(c_l, c_r)
do while (t < Time)
	
	t = t + dt
	
	do i = 1, N

		! call roeScheme(gamma, Cv, rho(i-1), u(i-1), p(i-1), rho(i), u(i), p(i), u_dot, Flux_minus)
		! call roeScheme(gamma, Cv, rho(i), u(i), p(i), rho(i+1), u(i+1), p(i+1), u_dot, Flux_plus)
		call hllScheme(gamma, Cv, rho(i-1), u(i-1), p(i-1), rho(i), u(i), p(i), u_dot, Flux_minus)
		call hllScheme(gamma, Cv, rho(i), u(i), p(i), rho(i+1), u(i+1), p(i+1), u_dot, Flux_plus)

		c_l = soundSpeed(gamma, rho(i-1), p(i-1), u(i-1))
		c_r = soundSpeed(gamma, rho(i+1), p(i+1), u(i+1))
		dt = CFL * dx  / max(c_l, c_r)
		
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
	u(0) = - u(1)
	p(0) = p(1)

	rho(N+1) = rho(N)
	u(N+1) = - u(N)
	p(N+1) = p(N)
	! end set boundary conditions
	
	!write(*,*) t

end do

call writeField(filename, N, x, rho, u, p)

end subroutine

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