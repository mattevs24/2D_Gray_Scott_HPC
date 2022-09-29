subroutine Mat_Mult(m,alpha,u,b)

  use header

  implicit none

  include 'mpif.h'

  !variables (inc internal variables)
  !=================================================
  integer, intent(in) :: m          !input n to multiply by
  real(kind=rk), intent(in):: alpha ! Multiple of matrix
  type(Vector), intent(in) :: u     ! Input vector u
  type(Vector), intent(inout) :: b  ! Output vector b = b + alpha*A*u

  integer :: i,ierr,nprocs,myid, snd_rqst1, recv_rqst1
  integer :: snd_rqst2,recv_rqst2, stat(MPI_STATUS_SIZE), nloc
  !=================================================
  
  !call the number of processors and processor rank
  !=================================================
  call MPI_comm_rank(MPI_COMM_WORLD,myid,ierr)
  call MPI_comm_size(MPI_COMM_WORLD,nprocs,ierr)
  !=================================================

  !collect and send the halo on either side of the u stored, using the non blocking methods
  !=================================================
  nloc = u%iend-u%ibeg+1

  call MPI_Isend(u%xx(nloc-m+1),m,MPI_DOUBLE_PRECISION,modulo(myid+1,nprocs),0, &
                                    MPI_COMM_WORLD,snd_rqst1,ierr)
  call MPI_Isend(u%xx(1),m,MPI_DOUBLE_PRECISION,modulo(myid-1,nprocs),1, &
                                    MPI_COMM_WORLD,snd_rqst2,ierr)
  call MPI_Irecv(u%xx(nloc+1),m,MPI_DOUBLE_PRECISION,modulo(myid-1,nprocs), &
                                    0, MPI_COMM_WORLD,recv_rqst1,ierr)
  call MPI_Irecv(u%xx(nloc+m+1),m,MPI_DOUBLE_PRECISION,modulo(myid+1,nprocs), &
                                    1, MPI_COMM_WORLD,recv_rqst2,ierr)
  !=================================================

  !old implementation as described in the writeup

  !do i = 1,nloc
     !top of box case
  !   if (modulo(i,m).EQ.1) then
  !      Aurow = -4*u%xx(i) + u%xx(i+1) + u%xx(i+m-1)
     !bottom of box case
  !   else if (modulo(i,m).EQ.0) then
  !      Aurow = -4*u%xx(i) + u%xx(i-1) + u%xx(i-m+1)
     !else case
  !   else
  !      Aurow = -4*u%xx(i) + u%xx(i+1) + u%xx(i-1)
  !   end if
  !   b%xx(i) = b%xx(i) + alpha * Aurow * m**2
  !end do

  !call MPI_wait(snd_rqst1,stat,ierr)
  !call MPI_wait(recv_rqst1,stat,ierr)
  !call MPI_wait(snd_rqst2,stat,ierr)
  !call MPI_wait(recv_rqst2,stat,ierr)

  !do i = 1,nloc
     !top of box case
  !   if (i.LE.m) then
  !      b%xx(i) = b%xx(i) + alpha * (u%xx(i+m)+u%xx(nloc+i)) * m**2
     !bottom of box case
  !   else if (i.GE.nloc-m+1) then
  !      b%xx(i) = b%xx(i) + alpha * (u%xx(i-m) + u%xx(2*m+i)) * m**2
  !   else
  !      b%xx(i) = b%xx(i) + alpha * (u%xx(i+m) + u%xx(i-m)) * m**2
  !   end if
  !end do

  !new update, calculate the b%xx(i) explicitly rather than in steps to save time and memory refences wrt Aurow
  !===================================================
  do i=1+m,nloc-m
     if (modulo(i,m).EQ.1) then
        b%xx(i) = b%xx(i) + alpha * m**2 * (-4*u%xx(i) + u%xx(i+1) + u%xx(i+m-1) + u%xx(i+m) + u%xx(i-m))
     else if (modulo(i,m).EQ.0) then
        b%xx(i) = b%xx(i) + alpha * m**2 * (-4*u%xx(i) + u%xx(i-1) + u%xx(i-m+1) + u%xx(i+m) + u%xx(i-m))
     else
        b%xx(i) = b%xx(i) + alpha * m**2 * (-4*u%xx(i) + u%xx(i-1) + u%xx(i+1) + u%xx(i+m) + u%xx(i-m))
     end if
  end do
  !This first part does the b%xx(i) that are in the middle of the local part of u%xx and hence don't use the...
  !data that is in the halos
  !===================================================

  !wait for the halo data to arrive from the non blocking calls earlier
  !===================================================
  call MPI_wait(snd_rqst1,stat,ierr)
  call MPI_wait(recv_rqst1,stat,ierr)
  call MPI_wait(snd_rqst2,stat,ierr)
  call MPI_wait(recv_rqst2,stat,ierr)
  !===================================================

  !uses the lower halo to calculate these b%xx(i)'s
  !===================================================
  do i = 1,m
     if(modulo(i,m).EQ.1) then
        b%xx(i) = b%xx(i) + alpha * m**2 * (-4*u%xx(i) + u%xx(i+1) + u%xx(i+m-1) + u%xx(i+m) + u%xx(nloc+i))
     else if (modulo(i,m).EQ.0) then
        b%xx(i) = b%xx(i) + alpha * m**2 * (-4*u%xx(i) + u%xx(i-1) + u%xx(i-m+1) + u%xx(i+m) + u%xx(nloc+i))
     else
        b%xx(i) = b%xx(i) + alpha * m**2 * (-4*u%xx(i) + u%xx(i-1) + u%xx(i+1) + u%xx(i+m) + u%xx(nloc+i))
     end if
  end do
  !===================================================

  !uses the upper halo to calculate these b%xx(i)'s
  !===================================================
  do i = nloc-m+1,nloc
     if(modulo(i,m).EQ.1) then
        b%xx(i) = b%xx(i) + alpha * m**2 * (-4*u%xx(i) + u%xx(i+1) + u%xx(i+m-1) + u%xx(i-m) + u%xx(2*m+i))
     else if (modulo(i,m).EQ.0) then
        b%xx(i) = b%xx(i) + alpha * m**2 * (-4*u%xx(i) + u%xx(i-1) + u%xx(i-m+1) + u%xx(i-m) + u%xx(2*m+i))
     else
        b%xx(i) = b%xx(i) + alpha * m**2 * (-4*u%xx(i) + u%xx(i-1) + u%xx(i+1) + u%xx(i-m) + u%xx(2*m+i))
     end if
  end do
  !===================================================
  !throughout this subroutine the modulo function is used to determine where in each block the i is. If...
  !it is in the very top slot then modulo(i,m)=1, very bottom then modulo(i,m)=0 and these have then the...
  !corner elements of the matrix included and shown in the diagrams on the assignment sheets

end subroutine Mat_Mult
