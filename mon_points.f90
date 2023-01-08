! finding cell inside which monitoring point lie
subroutine findMamaCell(N, x, surface_x, n_mon_point, x_mon_point, i_mon_point)
	implicit none

	integer										:: N, k, n_mon_point
	integer										:: i
	integer, dimension(n_mon_point)				:: i_mon_point
	double precision, dimension(n_mon_point)	:: x_mon_point
	double precision, dimension(0:N+1)			:: x
	double precision, dimension(N+1)			:: surface_x

	do k = 1, n_mon_point

		do i = 1, N

			if (x_mon_point(k) >= surface_x(i) .and. x_mon_point(k) <= surface_x(i+1)) then

				i_mon_point(k) = i

			end if

		end do

	end do

end subroutine