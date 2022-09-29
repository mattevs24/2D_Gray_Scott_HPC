subroutine RHS(m,u,v,F,k,Du,Dv,ru,rv)
  !==========================================
  ! Implement Gray-Scott right hand side here
  !==========================================

  use header

  implicit none

  !variables as needed
  !==========================================
  integer, intent(in) :: m
  type(Vector), intent(in) :: u,v
  real(kind=rk), intent(in) :: F,k,Du,Dv
  type(Vector), intent(inout) :: ru, rv

  !internal variables

  integer :: nloc
  !==========================================
  
  !calculate nloc on each processor
  !==========================================
  nloc = u%iend - u%ibeg + 1
  !==========================================
  
  !calculate rv and ru using both fortran's inbuilt functionality and daxpy
  !==========================================
  rv%xx = u%xx(1:nloc) * v%xx(1:nloc) * v%xx(1:nloc)
  ru%xx = -rv%xx + F
  
  call daxpy(nloc, -F, u%xx(1), 1, ru%xx(1), 1)
  call daxpy(nloc, -(F+k), v%xx(1), 1, rv%xx(1), 1)
  !==========================================

  !perform matmult on ru so far and the delta multiplication for the final ru and rv
  !==========================================
  call Mat_Mult(m,Du,u,ru)
  call Mat_Mult(m,Dv,v,rv)
  !==========================================

end subroutine RHS
