! HLL scheme for calculate fluxes vector
subroutine hllScheme(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, u_dot, Flux)
implicit none

integer 							:: i
double precision					:: rho_tilda, u_tilda, H_tilda, c_tilda, p_tilda, T_tilda, E_tilda
double precision					:: gamma, Cv
double precision					:: rho_l, u_l, H_l, p_l, c_l														! left parameters
double precision					:: rho_r, u_r, H_r, p_r, c_r 														! right parameters
double precision					:: enthalpy, soundSpeed																! functions declaration
double precision					:: u_dot																			! edge velo
double precision					:: S_l, S_r
double precision, dimension(3)		:: w_l, w_r, w_star
double precision, dimension(3)		:: F_l, F_r, F_star, Flux

H_l = enthalpy(gamma, rho_l, p_l, u_l)
H_r = enthalpy(gamma, rho_r, p_r, u_r)
c_l = soundSpeed(gamma, rho_l, p_l, u_l)
c_r = soundSpeed(gamma, rho_r, p_r, u_r)

! calculation tilda variables
rho_tilda = sqrt(rho_l * rho_r)
u_tilda = (sqrt(rho_l) * u_l + sqrt(rho_r) * u_r) / (sqrt(rho_l) + sqrt(rho_r))
H_tilda = (sqrt(rho_l) * H_l + sqrt(rho_r) * H_r) / (sqrt(rho_l) + sqrt(rho_r))

c_tilda = sqrt((gamma - 1.0) * (H_tilda - 0.5 * u_tilda**2))
p_tilda = (gamma - 1.0) / gamma * rho_tilda * (H_tilda - 0.5 * u_tilda**2)
T_tilda = 1 / Cv * (H_tilda - p_tilda / rho_tilda)
E_tilda = H_tilda - p_tilda / rho_tilda

call calcVariables(rho_l, u_l, p_l, H_l, w_l)
call calcVariables(rho_r, u_r, p_r, H_r, w_r)

call calcFluxes(rho_l, u_l, p_l, H_l, F_l)
call calcFluxes(rho_r, u_r, p_r, H_r, F_r)

!-------------HLL-------------
!S_l = min(u_l - c_l, u_r - c_r)
!S_r = max(u_l + c_l, u_r + c_r)

!-------------HLL(E)-------------
S_l = min(u_l - c_l, u_tilda - c_tilda)
S_r = max(u_r + c_r, u_tilda + c_tilda)

w_star = (S_r * w_r - S_l * w_l + F_l - F_r)/ (S_r - S_l)
F_star = (S_r * F_l - S_l * F_r + S_l * S_r * (w_r - w_l)) / (S_r - S_l)

if (S_l > u_dot) then
	Flux = F_l - u_dot * w_l
else if (S_r < u_dot) then
	Flux = F_r - u_dot * w_r
else
	Flux = F_star - u_dot * w_star
end if

end subroutine