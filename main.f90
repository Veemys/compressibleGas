program test
implicit none

character(*), parameter 				:: inputFile = "Input"
character(30) 							:: inputFileData

integer, parameter 						:: io = 12 				! in/out flag
integer									:: set_case
integer 								:: i, j, k

double precision						:: gamma, Cv
double precision						:: rho_l, u_l, H_l, p_l
double precision						:: rho_r, u_r, H_r, p_r
double precision						:: u_dot
double precision						:: enthalpy
double precision, dimension(3)			:: Flux

call test1D()

! write(*,*) "Enter 1 if you want to set parameters from console or 2 if you want to read from file..."
! read(*,*) set_case

! select case (set_case)
	! case (1)
		! call setParametersFromConsole(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, u_dot)
	! case (2)
		! open(io, file = inputFile)
			! read(io,*) inputFileData
		! close(io)
		! call readParametersFromFile(io, inputFileData, gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, u_dot)
	! case default
		! write(*,*) "Error input!!!"
		! stop
! end select

! H_l = enthalpy(gamma, rho_l, p_l, u_l)
! H_r = enthalpy(gamma, rho_r, p_r, u_r)

! call roeScheme(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, u_dot, Flux)

! call outputResult(Flux)

! Flux = 0.0
! call hllScheme(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, u_dot, Flux)

! call outputResult(Flux)

end