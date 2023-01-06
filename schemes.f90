! Roe scheme for calculate fluxes vector
subroutine roeScheme(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, u_dot, Flux)
	implicit none

	integer 							:: i
	double precision					:: rho_tilda, u_tilda, H_tilda, c_tilda, p_tilda, T_tilda, E_tilda
	double precision					:: lambda_1, lambda_2, lambda_3
	double precision					:: lambda_u_dot_1, lambda_u_dot_2, lambda_u_dot_3
	double precision					:: delta_u, delta_p, delta_rho
	double precision					:: gamma, Cv
	double precision					:: rho_l, u_l, H_l, p_l 														! left parameters
	double precision					:: rho_r, u_r, H_r, p_r 														! right parameters
	double precision					:: enthalpy																		! functions declaration
	double precision					:: u_dot																		! edge velo
	double precision, dimension(3)		:: v_temp1(3), v_temp2, v_temp3, deltaF_1, deltaF_2, deltaF_3
	double precision, dimension(3)		:: D, w_l, w_r, F_l, F_r, Flux
	double precision, dimension(3,3)	:: A

	! intent (in) gamma, rho_l, u_l, H_l, rho_r, u_r, H_r
	! intent (out) Flux

	! calculation enthalpy
	H_l = enthalpy(gamma, rho_l, p_l, u_l)
	H_r = enthalpy(gamma, rho_r, p_r, u_r)

	! calculation tilda variables
	rho_tilda = sqrt(rho_l * rho_r)
	u_tilda = (sqrt(rho_l) * u_l + sqrt(rho_r) * u_r) / (sqrt(rho_l) + sqrt(rho_r))
	H_tilda = (sqrt(rho_l) * H_l + sqrt(rho_r) * H_r) / (sqrt(rho_l) + sqrt(rho_r))

	c_tilda = sqrt((gamma - 1.0) * (H_tilda - 0.5 * u_tilda**2))
	p_tilda = (gamma - 1.0) / gamma * rho_tilda * (H_tilda - 0.5 * u_tilda**2)
	T_tilda = 1 / Cv * (H_tilda - p_tilda / rho_tilda)
	E_tilda = H_tilda - p_tilda / rho_tilda

	! calculation eigenvalues maxrix A
	lambda_1 = u_tilda - c_tilda
	lambda_2 = u_tilda
	lambda_3 = u_tilda + c_tilda

	lambda_u_dot_1 = max(u_dot, lambda_1) - min(u_dot, lambda_1)
	lambda_u_dot_2 = max(u_dot, lambda_2) - min(u_dot, lambda_2)
	lambda_u_dot_3 = max(u_dot, lambda_3) - min(u_dot, lambda_3)

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
	deltaF_1 = abs(lambda_u_dot_1) * ((delta_p - rho_tilda * c_tilda * delta_u) / (2.0 * c_tilda**2)) * v_temp1
	deltaF_2 = abs(lambda_u_dot_2) * (delta_rho - delta_p / c_tilda**2) * v_temp2
	deltaF_3 = abs(lambda_u_dot_3) * ((delta_p + rho_tilda * c_tilda * delta_u) / (2.0 * c_tilda**2)) * v_temp3

	! variables and flux'es calculate
	call calcVariables(rho_l, u_l, p_l, H_l, w_l)
	call calcVariables(rho_r, u_r, p_r, H_r, w_r)
	call calcFluxes(rho_l, u_l, p_l, H_l, F_l)
	call calcFluxes(rho_r, u_r, p_r, H_r, F_r)
	D = deltaF_1 + deltaF_2 + deltaF_3
	Flux = 0.5 * (F_l + F_r) - 0.5 * D - u_dot * 0.5 * (w_l + w_r)

end subroutine

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

! HLLC scheme for calculate fluxes vector
subroutine hllcScheme(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, Flux)
	implicit none

	integer 							:: i
	double precision					:: rho_tilda, u_tilda, H_tilda, c_tilda, p_tilda, T_tilda, E_tilda
	double precision					:: gamma, Cv
	double precision					:: rho_l, u_l, H_l, p_l, c_l														! left parameters
	double precision					:: rho_r, u_r, H_r, p_r, c_r 														! right parameters
	double precision					:: enthalpy, soundSpeed																! functions declaration
	double precision					:: S_l, S_r, S_star, p_star_l, p_star_r
	double precision, dimension(3)		:: D
	double precision, dimension(3)		:: w_l, w_r, w_star_l, w_star_r
	double precision, dimension(3)		:: F_l, F_r, F_star_l, F_star_r, Flux

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

	S_l = min(u_l - c_l, u_tilda - c_tilda)
	S_r = max(u_r + c_r, u_tilda + c_tilda)

	S_star = ((p_r - p_l) + rho_l * u_l * (S_l - u_l) - rho_r * u_r * (S_r - u_r)) / (rho_l * (S_l - u_l) - rho_r * (S_r - u_r))

	p_star_l = p_l + rho_l * (S_l - u_l) * (S_star - u_l)
	p_star_r = p_r + rho_r * (S_r - u_r) * (S_star - u_r)

	D(1) = 0
	D(2) = 1
	D(3) = S_star

	w_star_l = (S_l * w_l - F_l + p_star_l * D)/ (S_l - S_star)
	w_star_r = (S_r * w_r - F_r + p_star_r * D)/ (S_r - S_star)

	F_star_l = (S_star * (S_l * w_l - F_l) + S_l * p_star_l * D) / (S_l - S_star)
	F_star_r = (S_star * (S_r * w_r - F_r) + S_r * p_star_r * D) / (S_r - S_star)

	if (S_l > 0.0) then
		Flux = F_l
	else if (S_L <= 0.0 .and. S_star > 0.0) then
		Flux = F_star_l
	else if (S_star <= 0.0 .and. S_r > 0.0) then
		Flux = F_star_r
	else
		Flux = F_r
	end if

end subroutine

! Van Leer scheme for calculate fluxes vector
subroutine vanLeerScheme(gamma, Cv, rho_l, u_l, p_l, rho_r, u_r, p_r, Flux)
	implicit none

	integer 							:: i
	double precision					:: lambda_1, lambda_2, lambda_3
	double precision					:: gamma, Cv
	double precision					:: rho_l, u_l, H_l, p_l, c_l, M_l		! left parameters
	double precision					:: rho_r, u_r, H_r, p_r, c_r, M_r 		! right parameters
	double precision					:: coeff_l, coeff_r
	double precision					:: enthalpy
	double precision, dimension(3)		:: F_l, F_r, Flux

	! calculation enthalpy
	H_l = enthalpy(gamma, rho_l, p_l, u_l)
	H_r = enthalpy(gamma, rho_r, p_r, u_r)

	! calculation of sound speed and Mach number
	c_l = sqrt(gamma * p_l / rho_l)
	c_r = sqrt(gamma * p_r / rho_r)
	M_l = u_l / c_l
	M_r = u_r / c_r

	! calculation eigenvalues maxrix A
	lambda_1 = 1.0
	lambda_2 = 1.0
	lambda_3 = 1.0

	! left and right fluxes calculation
	coeff_l = rho_l * c_l / 4.0 * (M_l + 1.0)**2
	coeff_r = - rho_r * c_r / 4.0 * (M_r - 1.0)**2
	F_l(1) = coeff_l * 1.0
	F_l(2) = coeff_l * 2.0 * c_l / gamma * (1.0 + (gamma - 1.0) / 2.0 * M_l)
	F_l(3) = coeff_l * 2.0 * c_l**2 / (gamma**2 - 1.0) * (1.0 + (gamma - 1.0) / 2.0 * M_l)**2

	F_r(1) = coeff_r * 1.0
	F_r(2) = coeff_r * 2.0 * c_r / gamma * (-1.0 + (gamma - 1.0) / 2.0 * M_r)
	F_r(3) = coeff_r * 2.0 * c_r**2 / (gamma**2 - 1.0) * (1.0 - (gamma - 1.0) / 2.0 * M_r)**2

	! flux'es calculate
	Flux = F_l + F_r

end subroutine


! TVD scheme
subroutine TVD_MUSCL(i, N, rho, u, p, &
					 rho_ll, u_ll, p_ll, rho_lr, u_lr, p_lr, &
					 rho_rl, u_rl, p_rl, rho_rr, u_rr, p_rr)
	implicit none

	integer									:: i, N
	double precision						:: rho_ll, u_ll, p_ll, rho_lr, u_lr, p_lr
	double precision						:: rho_rl, u_rl, p_rl, rho_rr, u_rr, p_rr
	double precision, dimension(0:N + 1)	:: rho, u, p
	double precision						:: limiter

	if (i == 1) then

		rho_ll = rho(0)
		u_ll = u(0)
		p_ll = p(0)

		rho_lr = rho(i) - 0.5 * limiter(rho(i-1), rho(i), rho(i+1)) * (rho(i) - rho(i - 1))
		u_lr = u(i) - 0.5 * limiter(u(i-1), u(i), u(i+1)) * (u(i) - u(i - 1))
		p_lr = p(i) - 0.5 * limiter(p(i-1), p(i), p(i+1)) * (p(i) - p(i - 1))

		rho_rl = rho(i) + 0.5 * limiter(rho(i-1), rho(i), rho(i+1)) * (rho(i) - rho(i-1))
		u_rl = u(i) + 0.5 * limiter(u(i-1), u(i), u(i+1)) * (u(i) - u(i-1))
		p_rl = p(i) + 0.5 * limiter(p(i-1), p(i), p(i+1)) * (p(i) - p(i-1))

		rho_rr = rho(i+1) - 0.5 * limiter(rho(i), rho(i+1), rho(i+2)) * (rho(i+1) - rho(i))
		u_rr = u(i+1) - 0.5 * limiter(u(i), u(i+1), u(i+2)) * (u(i+1) - u(i))
		p_rr = p(i+1) - 0.5 * limiter(p(i), p(i+1), p(i+2)) * (p(i+1) - p(i))

	else if (i == N) then

		rho_ll = rho(i-1) + 0.5 * limiter(rho(i-2), rho(i-1), rho(i)) * (rho(i-1) - rho(i-2))
		u_ll = u(i-1) + 0.5 * limiter(u(i-2), u(i-1), u(i)) * (u(i-1) - u(i-2))
		p_ll = p(i-1) + 0.5 * limiter(p(i-2), p(i-1), p(i)) * (p(i-1) - p(i-2))

		rho_lr = rho(i) - 0.5 * limiter(rho(i-1), rho(i), rho(i+1)) * (rho(i) - rho(i - 1))
		u_lr = u(i) - 0.5 * limiter(u(i-1), u(i), u(i+1)) * (u(i) - u(i - 1))
		p_lr = p(i) - 0.5 * limiter(p(i-1), p(i), p(i+1)) * (p(i) - p(i - 1))

		rho_rl = rho(i) + 0.5 * limiter(rho(i-1), rho(i), rho(i+1)) * (rho(i) - rho(i-1))
		u_rl = u(i) + 0.5 * limiter(u(i-1), u(i), u(i+1)) * (u(i) - u(i-1))
		p_rl = p(i) + 0.5 * limiter(p(i-1), p(i), p(i+1)) * (p(i) - p(i-1))

		rho_rr = rho(n+1)
		u_rr = u(n+1)
		p_rr = p(n+1)

	else

		rho_ll = rho(i-1) + 0.5 * limiter(rho(i-2), rho(i-1), rho(i)) * (rho(i-1) - rho(i-2))
		u_ll = u(i-1) + 0.5 * limiter(u(i-2), u(i-1), u(i)) * (u(i-1) - u(i-2))
		p_ll = p(i-1) + 0.5 * limiter(p(i-2), p(i-1), p(i)) * (p(i-1) - p(i-2))

		rho_lr = rho(i) - 0.5 * limiter(rho(i-1), rho(i), rho(i+1)) * (rho(i) - rho(i-1))
		u_lr = u(i) - 0.5 * limiter(u(i-1), u(i), u(i+1)) * (u(i) - u(i - 1))
		p_lr = p(i) - 0.5 * limiter(p(i-1), p(i), p(i+1)) * (p(i) - p(i - 1))

		rho_rl = rho(i) + 0.5 * limiter(rho(i-1), rho(i), rho(i+1)) * (rho(i) - rho(i-1))
		u_rl = u(i) + 0.5 * limiter(u(i-1), u(i), u(i+1)) * (u(i) - u(i-1))
		p_rl = p(i) + 0.5 * limiter(p(i-1), p(i), p(i+1)) * (p(i) - p(i-1))

		rho_rr = rho(i+1) - 0.5 * limiter(rho(i), rho(i+1), rho(i+2)) * (rho(i+1) - rho(i))
		u_rr = u(i+1) - 0.5 * limiter(u(i), u(i+1), u(i+2)) * (u(i+1) - u(i))
		p_rr = p(i+1) - 0.5 * limiter(p(i), p(i+1), p(i+2)) * (p(i+1) - p(i))

	end if

end subroutine

double precision function limiter(u_prev, u, u_next)
	implicit none
	
	double precision				:: r
	double precision				:: u_prev, u, u_next

	r = (u_next - u) / (u - u_prev + 1e-16)
	if (r <= 0) then

		limiter = 0
		return

	end if

	! limiter = min(1.0, r)									! minmod
	! limiter = (r**2 + r) / (r**2 + 1)						! van Albada
	! limiter = 2.0 * r / (r + 1.0)							! van Leer
	limiter = max(min(2.0 * r, 1.0), min(r, 2.0))			! superbee

end function

! вадим
! subroutine TVD_MUSCL(R, P, U, n, R_f_tvd, P_f_tvd, U_f_tvd)
! 	implicit none

! 	integer n, i
! 	Real(8) R(0:n+1), P(0:n+1), U(0:n+1)
! 	Real(8) R_f_tvd(1:n+1,1:2), P_f_tvd(1:n+1,1:2), U_f_tvd(1:n+1,1:2)
! 	Real(8) r_r, r_p, r_u
! 	Real(8) limiter

! 	R_f_tvd(1,1) = R(0) ! пока определим значнеие слева на первую границу, как просто значение из заграничной ячейки
! 	P_f_tvd(1,1) = P(0)
! 	U_f_tvd(1,1) = U(0)
! 	R_f_tvd(n+1,2) = R(n+1) ! аналогично для последней грани
! 	P_f_tvd(n+1,2) = P(n+1)
! 	U_f_tvd(n+1,2) = U(n+1)

! 	do i = 1 , n
	
! 		R_f_tvd(i,2) = R(i) - 0.5*limiter(R(i+1)-R(i), R(i)-R(i-1) )*(R(i)-R(i-1))
! 		R_f_tvd(i+1,1) = R(i) + 0.5*limiter(R(i+1)-R(i), R(i)-R(i-1) )*(R(i)-R(i-1))

! 		P_f_tvd(i,2) = P(i) - 0.5*limiter(P(i+1)-P(i) , P(i)-P(i-1) )*(P(i)-P(i-1))
! 		P_f_tvd(i+1,1) = P(i) + 0.5*limiter(P(i+1)-P(i) , P(i)-P(i-1) )*(P(i)-P(i-1))

! 		U_f_tvd(i,2) = U(i) - 0.5*limiter(U(i+1)-U(i), U(i)-U(i-1))*(U(i)-U(i-1))
! 		U_f_tvd(i+1,1) = U(i) + 0.5*limiter( U(i+1)-U(i), U(i)-U(i-1) )*(U(i)-U(i-1))
		
! 	end do

! end subroutine

! function limiter(r1, r2 )
! 	implicit none
! 	Real(8) limiter, r , r1 , r2

! 	if (r2 < 1.e-12 ) then
	
! 		limiter=0.
! 		return
		
! 	else

! 		r=r1/r2
! 		if (r<1e-12) then
! 			limiter=0.
! 		else
! 			limiter=(r**2 + r)/(r**2 + 1.0) ! симметричный
! 			!limiter = (2.*r)/(r+1.) ! тоже симметричный
! 		end if
		
! 	end if

! end function