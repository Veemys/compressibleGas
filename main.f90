program test
	implicit none

	character(*), parameter 				:: inputFile = "Input"
	character(30) 							:: inputFileData

	integer, parameter 						:: io = 12 				! in/out flag
	integer									:: set_case
	integer 								:: i, j, k

	integer									:: N

	double precision						:: gamma, Cv
	double precision						:: L, x_0, Time, CFL
	double precision						:: rho_l, u_l, p_l, H_l
	double precision						:: rho_r, u_r, p_r, H_r
	double precision						:: rho, u, p
	double precision						:: u_dot
	double precision						:: enthalpy
	double precision, dimension(3)			:: Flux
	
	! **********************************
	! *			 расчет трубы		   *
	! **********************************
	
	open(io, file = inputFileData)
	
	read(io,*) L				 	! Domain length
	read(io,*) x_0					! Initial discontinuity position
	read(io,*) N					! Number of computing cells
	read(io,*) gamma				! Ratio of specific heats
	read(io,*) Cv					! Cv
	read(io,*) Time					! Output time
	read(io,*) CFL					! Courant number
	read(io,*) rho_l				! Initial density on left state
	read(io,*) u_l					! Initial velocity on left state
	read(io,*) p_l					! Initial pressure on left state
	read(io,*) rho_r				! Initial density on right state
	read(io,*) u_r					! Initial velocity on right state
	read(io,*) p_r					! Initial pressure on right

	close(io)
	
	write(*,*) "L = ", L
	write(*,*) "x0 = ", x_0
	write(*,*) "N = ", N
	write(*,*) "gamma = ", gamma
	write(*,*) "Cp = ", Cv
	write(*,*) "Time = ", Time
	write(*,*) "CFL = ", CFL
	write(*,*) "rho_l/rho_r = ", rho_l, "/", rho_r
	write(*,*) "u_l/u_r = ", u_l, "/", u_r
	write(*,*) "p_l/p_r = ", p_l, "/", p_r

	! call test1D(L, x_0, N, gamma, Cv, Time, CFL, rho_l, u_l, p_l, rho_r, u_r, p_r)

	call FVM(L, x_0, N, gamma, Cv, Time, CFL, rho_l, u_l, p_l, rho_r, u_r, p_r)

	! **********************************
	! *	тест расчета потоков на гранях *
	! **********************************
	
	! open(io, file = inputFile)
	! read(io,*) inputFileData
	! close(io)

	! call readParametersFromFile(io, inputFileData, gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, u_dot)

	! H_l = enthalpy(gamma, rho_l, p_l, u_l)
	! H_r = enthalpy(gamma, rho_r, p_r, u_r)

	! call roeScheme(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, u_dot, Flux)

	! call outputResult(Flux)

	! Flux = 0.0
	! call hllScheme(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, u_dot, Flux)

	! call outputResult(Flux)

end