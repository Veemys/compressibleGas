! Roe scheme for calculate fluxes vector
subroutine roeScheme(gamma, Cv, rho_l, u_l, p_l, H_l, rho_r, u_r, p_r, H_r, Flux)
implicit none

integer 							:: i
double precision					:: rho_tilda, u_tilda, H_tilda, c_tilda, p_tilda, T_tilda, E_tilda
double precision					:: lambda_1, lambda_2, lambda_3
double precision					:: delta_u, delta_p, delta_rho
double precision					:: gamma, Cv
double precision					:: rho_l, u_l, H_l, p_l 														! left parameter
double precision					:: rho_r, u_r, H_r, p_r 														! right parameter
double precision, dimension(3)		:: v_temp1(3), v_temp2, v_temp3, deltaF_1, deltaF_2, deltaF_3
double precision, dimension(3)		:: D, F_l, F_r, Flux
double precision, dimension(3,3)	:: A

! intent (in) gamma, rho_l, u_l, H_l, rho_r, u_r, H_r
! intent (out) Flux

! calculation tilda variables
rho_tilda = sqrt(rho_l * rho_r)
u_tilda = (sqrt(rho_l) * u_l + sqrt(rho_r) * u_r) / (sqrt(rho_l) + sqrt(rho_r))
H_tilda = (sqrt(rho_l) * H_l + sqrt(rho_r) * H_r) / (sqrt(rho_l) + sqrt(rho_r))

c_tilda = sqrt((gamma - 1.0) * (H_tilda - 0.5 * u_tilda**2))
p_tilda = (gamma - 1.0) / gamma * rho_tilda * (H_tilda - 0.5 * u_tilda**2)
T_tilda = 1 / Cv * (H_tilda - p_tilda / rho_tilda)
E_tilda = H_tilda - p_tilda / rho_tilda

! calculation eigenvalues maxrix A
lambda_1 = u_tilda
lambda_2 = u_tilda + c_tilda
lambda_3 = u_tilda - c_tilda

! calc matrix A with Roe variables
call calcAMatrix(gamma, u_tilda, rho_tilda, E_tilda, A)

! delta u, p, rho calculation
delta_u = u_r - u_l
delta_p = p_r - p_l
delta_rho = rho_r - rho_l

! temporary vectors for calculation deltaF
v_temp1(1) = 1
v_temp1(2) = u_tilda - c_tilda
v_temp1(3) = H_tilda - c_tilda * u_tilda

v_temp2(1) = 1
v_temp2(2) = u_tilda
v_temp2(3) = u_tilda**2 / 2.0

v_temp3(1) = 1
v_temp3(2) = u_tilda + c_tilda
v_temp3(3) = H_tilda + c_tilda * u_tilda

! deltaF calculation
deltaF_1 = abs(u_tilda - c_tilda) * ((delta_p - rho_tilda * c_tilda * delta_u) / (2.0 * c_tilda**2)) * v_temp1
deltaF_2 = abs(u_tilda) * (delta_rho - delta_p / c_tilda**2) * v_temp2
deltaF_3 = abs(u_tilda + c_tilda) * ((delta_p + rho_tilda * c_tilda * delta_u) / (2.0 * c_tilda**2)) * v_temp3

! flux'es calculate
call calcFluxes(rho_l, u_l, p_l, H_l, F_l)
call calcFluxes(rho_r, u_r, p_r, H_r, F_r)
D = deltaF_1 + deltaF_2 + deltaF_3
Flux = 0.5 * (F_l + F_r) - 0.5 * D

end subroutine

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