! wall
subroutine BC_wall(rho_in, u_in, p_in, rho_out, u_out, p_out)
	implicit none

	double precision 	:: rho_in, u_in, p_in, rho_out, u_out, p_out

	rho_out = rho_in
	u_out = - u_in
	p_out = p_in

end subroutine

! transmissive
subroutine BC_transmissive(rho_in, u_in, p_in, rho_out, u_out, p_out)
	implicit none
	
	double precision	:: rho_in, u_in, p_in, rho_out, u_out, p_out
	
	rho_out = rho_in
	u_out = u_in
	p_out = p_in

end subroutine

! piston
subroutine BC_piston(rho_in, u_in, p_in, rho_out, u_out, p_out, u_p)
	implicit none
	
	double precision	:: rho_in, u_in, p_in, rho_out, u_out, p_out, u_p 	! u_p - sound of piston
	
	rho_out = rho_in
	u_out = 2.0 * u_p - u_in
	p_out = p_in

end subroutine

! open boundary
subroutine BC_open(gamma, rho_in, u_in, p_in, rho_out, u_out, p_out, norm, rho_0, P_0)
	implicit none

	double precision	:: rho_in, u_in, p_in, rho_out, u_out, p_out
	double precision	:: norm, rho_0, P_0
	double precision	:: gamma, c, u_n
	double precision	:: soundSpeed

	u_n = u_in * norm
	c = soundSpeed(gamma, rho_in, p_in, u_in)
	
	! inlet
	if (u_n > 0) then
		
		! supersonic
		if (abs(u_n) > c) then
			
			rho_out = rho_0 * (1.0 + (gamma - 1.0) / 2.0 * (u_in / c)**2)**(- 1.0 / (gamma - 1.0))
			u_out = c
			p_out = P_0 * (1.0 + (gamma - 1.0) / 2.0 * (u_in / c)**2)**(- gamma / (gamma - 1.0))
			
		! subsonic
		else
		
			rho_out = rho_0 * (1.0 + (gamma - 1.0) / 2.0 * (u_in / c)**2)**(- 1.0 / (gamma - 1.0))
			u_out = u_in
			p_out = P_0 * (1.0 + (gamma - 1.0) / 2.0 * (u_in / c)**2)**(- gamma / (gamma - 1.0))
			
		end if
		
	! outlet
	else
	
		! supersonic
		if (abs(u_n) > c) then
		
			rho_out = rho_in
			u_out = u_in
			p_out = p_in
			
		! subsonic
		else
		
			rho_out = rho_in
			u_out = u_in
			p_out = P_0
			
		end if
		
	end if

end subroutine