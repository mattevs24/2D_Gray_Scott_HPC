subroutine initial(m, u, v)
  !========================================
  ! A subroutine for assembling 2-Gaussians
  ! initial condition
  !========================================
  use header
  implicit none

  ! Variables as per assignment
  integer, intent(in) :: m
  type(Vector), intent(inout) :: u,v

  ! Loop counters: total, horizontal, vertical
  integer ::       irow,  i,          j, k
  ! mesh step, horizontal, vertical coordinates
  real(kind=rk) :: h, x, y, vi
  ! 2*variance of Gaussians

  real(kind=rk), parameter :: v2 = 0.02_rk
  
  !added the additional variable vi (v inverse) here to save on the divisions later in the initial code
  h = 1.0_rk/m
  vi = 1.0_rk/v2
      
  do irow=u%ibeg,u%iend
     ! Calculate Cartesian index splitting

     !needed to rescale and so used k to get the actual start and end points of the vector input
     j = (irow-1)/m + 1
     i = irow - (j-1)*m
     y = j*h
     x = i*h
     k = irow-u%ibeg+1

     u%xx(k) = 1.0_rk - 0.5_rk*exp(-vi*(x-0.5_rk)**2 -vi*(y-0.5_rk)**2) &
                         - 0.5_rk*exp(-vi*(x-0.4_rk)**2 -vi*(y-0.6_rk)**2)
     v%xx(k) = 0.25_rk*exp(-vi*(x-0.4_rk)**2 -vi*(y-0.4_rk)**2)  &
                + 0.25_rk*exp(-vi*(x-0.5_rk)**2 -vi*(y-0.6_rk)**2)     
  end do

end subroutine initial
     
