program grayscott
  ! ===================================
  ! Main program for solving
  ! 2-dimensional Gray-Scott equations
  ! (sequential version)
  ! ===================================
  
  use header

  implicit none

  include "mpif.h"

  ! Variables as per assignment
  !=======================================================
  integer :: m, ierr, myid, nprocs, ibeg, iend, nrows
  type(Vector) :: u, v 
  real(kind=rk) :: tau, T, F, k, Du, Dv

  ! 2D problem size
  integer :: n
  ! Ticks for CPU timing
  real(kind=rk) :: t_start, t_finish
  type(vector) :: utotal, vtotal
  !=======================================================

  !initialise the MPI and call in the comm rank and size
  !=======================================================
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
  !=======================================================
    
  !initialise the data from the input.dat file
  !=======================================================
  if (myid == 0) then
     open(unit=2,file="input.dat")
     read(2,*) m
     read(2,*) tau
     read(2,*) T
     read(2,*) F
     read(2,*) k
     read(2,*) Du
     read(2,*) Dv
     close(2)
     write(*, '(A,I5)')   'm  = ', m
     write(*, '(A,F7.3)') 'tau= ', tau
     write(*, '(A,F8.1)') 'T  = ', T
     write(*, '(A,F9.5)') 'F  = ', F
     write(*, '(A,F9.5)') 'k  = ', k
     write(*, '(A,ES11.3)')'Du = ', Du
     write(*, '(A,ES11.3)')'Dv = ', Dv
  end if
  !=======================================================
  
  !broadcast this data to the other processors
  !=======================================================
  call MPI_Bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tau,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(T,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(F,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(k,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Du,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Dv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !=======================================================

  !Calculate n and set all the vectors n
  !=======================================================
  n = m*m !in all processors

  !processor specific ibeg,iend and nrows
  ibeg = nint(myid*n/real(nprocs)) + 1
  iend = nint((myid+1)*n/real(nprocs))
  nrows = iend - ibeg + 1

  ! Initial guess
  u%n = n
  ! Sizes are important for correct allocation in initial
  v%n = n
  !=======================================================

  !allocate and occupy u and v with ibeg, iend and xx
  !=======================================================
  allocate(u%xx(nrows+2*m),v%xx(nrows+2*m))
   
  u%ibeg = ibeg
  u%iend = iend

  v%ibeg = ibeg
  v%iend = iend

  call initial(m,u,v)
  !=======================================================

  ! Run the Heun method
  !=======================================================
  call cpu_time(t_start)
  call timestepping(u,v, tau,T,F,k,Du,Dv, m)
  call cpu_time(t_finish)
  !=======================================================

  !write the time and the results for that one value in u
  !=======================================================
  if (myid == nprocs-1) then
     write(*, '(A,F17.12)') 'u_{n/2,n}     = ', u%xx(nrows-m/2)
     write(*, '(A,F8.4)')  'total cpu time = ', t_finish - t_start
  end if
  !=======================================================

  !on the first processor only, allocate the u and v total vectors
  !=======================================================
  if (myid==0) then
     allocate(utotal%xx(n),vtotal%xx(n))
     utotal%n = n
     vtotal%n = n
  end if
  !=======================================================

  !Gather all the parts of u and v together in u/vtotal
  !=======================================================
  call MPI_Gather(u%xx(1),nrows,MPI_DOUBLE_PRECISION,utotal%xx(nrows*myid+1),nrows,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Gather(v%xx(1),nrows,MPI_DOUBLE_PRECISION,vtotal%xx(nrows*myid+1),nrows,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !=======================================================

  !complete the save_fields call and save u and v to solution.dat and then deallocate utotal%xx and vtotal%xx
  !=======================================================
  if (myid==0) then
    call save_fields(utotal,vtotal,'solution.dat')
    deallocate(utotal%xx,vtotal%xx)
  end if
  !=======================================================

  !deallocate the local u and v parts on all the processors
  !=======================================================
  deallocate(u%xx, v%xx)
  !=======================================================

  !finalise the MPI
  !=======================================================
  call MPI_Finalize(ierr)
  !=======================================================

end program grayscott
