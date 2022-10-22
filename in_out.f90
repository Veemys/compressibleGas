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