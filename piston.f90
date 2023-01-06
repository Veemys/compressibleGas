! calc piston's velocity
double precision function veloPiston(amplitude, frequency, time)
	implicit none
	
	double precision	:: amplitude, frequency, time
	
	veloPiston = amplitude * frequency * sin(frequency * time)
	
end function

! calc new piston's position
double precision function xPiston(amplitude, frequency, time)
	implicit none
	
	double precision	::  amplitude, frequency, time
	
	xPiston = - amplitude * cos(frequency * time)
	! xPiston = x_old + u_piston * dt
	
end function