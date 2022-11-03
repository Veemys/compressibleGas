! calculation of total enthalpy
double precision function enthalpy(gamma, rho, p, u)
implicit none

double precision 	:: gamma, rho, p, u

enthalpy = p / rho * (gamma / (gamma - 1.0)) + u**2 / 2.0

end function

! calculation of speed of sound
double precision function soundSpeed(gamma, rho, p, u)
implicit none

double precision	:: gamma, rho, p, u

soundSpeed = sqrt(gamma * p / rho)

end function


! calculate A matrix
subroutine calcAMatrix(gamma, u, rho, E, A)
implicit none

double precision					:: gamma, u, rho, E
double precision, dimension(3,3)	:: A
intent (in)  gamma, u, rho, E
intent (out) A

A(1,1) = 0
A(1,2) = 1
A(1,3) = 0

A(2,1) = (gamma - 3.0) * u**2 / 2.0
A(2,2) = (gamma - 3.0) * u
A(2,3) = (gamma - 1.0)

A(3,1) = (gamma - 1.0) * u**3 - gamma * u * E / rho
A(3,2) = gamma * E / rho - 3.0 * (gamma - 1.0) * u**2 / 2.0
A(3,3) = gamma * u

end subroutine

! calculation fluxes vector
subroutine calcFluxes(rho, u, p, H, F)
implicit none

double precision					:: rho, u, p, H
double precision, dimension(3)		:: F
intent (in)  rho, u, p, H
intent (out) F

F(1) = rho * u
F(2) = rho * u**2 + p
F(3) = rho * H * u

end subroutine

! calculation variables vector
subroutine calcVariables(rho, u, p, H, w)
implicit none

double precision					:: rho, u, p, H
double precision, dimension(3)		:: w
intent (in)  rho, u, p, H
intent (out) w

w(1) = rho
w(2) = rho * u
w(3) = rho * (H - p / rho)

end subroutine

! calculation variables from variables vector
subroutine calcVariablesFromVector(gamma, rho, u, p, w)
implicit none

double precision					:: gamma, rho, u, p
double precision, dimension(3)		:: w
intent (in)  w
intent (out) rho, u, p

rho = w(1)
u = w(2)/ rho
p = (rho * u**2 / 2.0 - w(3)) / (1.0 - gamma / (gamma - 1.0))

end subroutine