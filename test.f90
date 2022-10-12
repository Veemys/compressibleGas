program test
implicit none

integer 								:: io ! in/out flag
integer 								:: i, j, k
double precision						:: gamma, Cv
double precision						:: rho_l, u_l, H_l, p_l
double precision						:: rho_r, u_r, H_r, p_r
double precision						:: enthalpy
double precision, dimension(3)			:: Flux

call setParametersFromConsole(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r)

H_l = enthalpy(gamma, rho_l, p_l, u_l)
H_r = enthalpy(gamma, rho_r, p_r, u_r)

call roeScheme(gamma, Cv, rho_l, u_l, p_l, H_l, rho_r, u_r, p_r, H_r, Flux)

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

character			:: filename
integer				:: io
double precision 	:: gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r

open(io, file = filename)
read(io,*) gamma
read(io,*) Cv
read(io,*) rho_l, rho_r
read(io,*) u_l, u_r
read(io,*) p_l, p_r
close(io)

end subroutine

subroutine outputResult(Flux)
implicit none

double precision, dimension(3) 	:: Flux

write(*,*) "Fluxes = ", Flux(1), Flux(2), Flux(3)

end subroutine

! calculation of total enthalpy
double precision function enthalpy(gamma, rho, p, u)
implicit none

double precision 	:: gamma, rho, p, u

enthalpy = p / rho * (gamma / (gamma - 1.0)) + u**2 / 2.0

end function