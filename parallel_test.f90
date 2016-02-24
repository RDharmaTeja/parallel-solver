! test the parallel module
program test
use mpi
! use mpi.h if mpi module is not available
use parallel
implicit none
integer :: ierr, current_id,numprocs
real :: buffer
integer :: status(MPI_STATUS_SIZE)

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,current_id,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

print *, "current rank is ", current_id
call read_layout_file(current_id)

if(current_id == 0) then
buffer = 111.1
call MPI_SEND(buffer,1,MPI_REAL,1,1,MPI_COMM_WORLD, ierr)
call MPI_RECV(buffer,1,MPI_REAL,1,1,MPI_COMM_WORLD,status,ierr)
print *, "received from 1 ", buffer
else
call MPI_RECV(buffer,1,MPI_REAL,0,1,MPI_COMM_WORLD,status,ierr)
buffer = 111.122
call MPI_SEND(buffer,1,MPI_REAL,0,1,MPI_COMM_WORLD, ierr)
print *, "received from 1 ", buffer
end if

call MPI_FINALIZE(ierr)

end program test
