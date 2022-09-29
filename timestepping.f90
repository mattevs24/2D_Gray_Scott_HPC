subroutine timestepping(u,v,tau,T,F,k,Du,Dv,m)
  !=========================================
  ! Implement Heun method (Algorithm 1) here
  !=========================================
  
  use header

  implicit none

  !variables as per the assignment
  !=========================================
  type(Vector), intent(inout) :: u,v
  integer, intent(in) :: m
  real(kind=rk), intent(in) :: tau, T, F, k, Du, Dv

  !internal use variables here
  real(kind=rk) :: tauhalf !speeds up divisions
  type(Vector) :: ru, ruhat
  type(Vector) :: rv, rvhat
  integer :: i, nloc
  !=========================================

  !calculating the tauhalf and nloc variables to save time and memory references later
  !=========================================
  tauhalf = tau/2
  nloc = u%iend - u%ibeg + 1
  !=========================================

  !allocate the ru/v, ruhat,rvhat vectors
  !=========================================
  allocate(ru%xx(nloc),rv%xx(nloc),ruhat%xx(nloc),rvhat%xx(nloc))
  !=========================================

  !perform algorithm 1 (heun method here)
  !=========================================
  do i = 0,(nint(T/tau)-1)

     !Step 1
     call RHS(m,u,v,F,k,Du,Dv,ru,rv)
    
     !Step 2
     call daxpy(nloc,tau,ru%xx(1),1,u%xx(1),1)
     call daxpy(nloc,tau,rv%xx(1),1,v%xx(1),1)

     !Step 3
     call RHS(m,u,v,F,k,Du,Dv,ruhat,rvhat)

     !Step 4
     call daxpy(nloc,-1.0_8,ru%xx(1),1,ruhat%xx(1),1)
     call daxpy(nloc,-1.0_8,rv%xx(1),1,rvhat%xx(1),1)

     call daxpy(nloc,tauhalf,ruhat%xx(1),1,u%xx(1),1)
     call daxpy(nloc,tauhalf,rvhat%xx(1),1,v%xx(1),1)
     
  end do
  !=========================================

  !deallocation
  !=========================================
  deallocate(ru%xx,rv%xx,ruhat%xx,rvhat%xx)
  !=========================================

end subroutine timestepping

