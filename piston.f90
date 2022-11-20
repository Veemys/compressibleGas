! calc piston's velocity
double precision function veloPiston(amplitude, frequency, time)
	implicit none
	
	double precision							:: amplitude, frequency, time
	
	veloPiston = amplitude * frequency * sin(frequency * time)
	
end function

! calc new piston's position
double precision function xPiston(x_old, u_piston, dt)
	implicit none
	
	double precision							::  x_old, u_piston, dt
	
	xPiston = x_old + u_piston * dt
	
end function