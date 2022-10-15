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
double precision						:: enthalpy
double precision, dimension(3)			:: Flux

write(*,*) "Enter 1 if you want to set parameters from console or 2 if you want to read from file..."
read(*,*) set_case

select case (set_case)
	case (1)
		call setParametersFromConsole(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r)
	case (2)
		open(io, file = inputFile)
			read(io,*) inputFileData
		close(io)
		call readParametersFromFile(io, inputFileData, gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r)
	case default
		write(*,*) "Error input!!!"
		stop
end select

H_l = enthalpy(gamma, rho_l, p_l, u_l)
H_r = enthalpy(gamma, rho_r, p_r, u_r)

call roeScheme(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, Flux)

call outputResult(Flux)

end

! setting parameters from console
subroutine setParametersFromConsole(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r)
implicit none

double precision 	:: gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r

write(*,*) "Enter a gamma: "
read(*,*) gamma
write(*,*) "Enter a thermal capacity Cv: "
read(*,*) Cv
write(*,*) "Enter a left and right density: "
read(*,*) rho_l, rho_r
write(*,*) "Enter a left and right velocity: "
read(*,*) u_l, u_r
write(*,*) "Enter a left and right pressure: "
read(*,*) p_l, p_r

end subroutine

! reading parameters from file
subroutine readParametersFromFile(io, filename, gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r)
implicit none

character(30)		:: filename
integer				:: io
double precision 	:: gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r
double precision	:: u_dot

open(io, file = filename)
read(io,*) gamma
read(io,*) Cv
read(io,*) p_l, p_r
read(io,*) rho_l, rho_r
read(io,*) u_l, u_r
read(io,*) u_dot
close(io)

end subroutine

! writing output to console
subroutine outputResult(Flux)
implicit none

double precision, dimension(3) 	:: Flux

write(*,*) "Fluxes = ", Flux

end subroutine

! calculation of total enthalpy
double precision function enthalpy(gamma, rho, p, u)
implicit none

double precision 	:: gamma, rho, p, u

enthalpy = p / rho * (gamma / (gamma - 1.0)) + u**2 / 2.0

end function