subroutine save_fields(u,v,filename)
  ! ===========================================================
  ! Save solutions to disk.
  ! Writes a text file which contains the values of the vectors
  ! u and v at the nodal points.
  ! ===========================================================

  use header
  implicit none
  type(vector), intent(in) :: u, v
  ! Name of file to write to
  character(len=*), intent(in) :: filename
  ! Loop variable
  integer :: i
  integer :: file_id

  file_id = 99 ! avoid conflicts if file of parameters is not closed

  !Write first component
  open(unit=file_id,file=trim(filename))
  do i=1, u%n
     write(file_id,'(E20.8e3," ")',advance='no') u%xx(i)
  end do
  write(file_id,'("")') ! delimiter
  !Write second component
  do i=1, v%n
     write(file_id,'(E20.8e3," ")',advance='no') v%xx(i)
  end do
  write(file_id,'("")')

  close(file_id)

end subroutine save_fields
