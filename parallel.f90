module parallel
  ! contains routines to execute parallelly
  use mpi
  use global, only: CONFIG_FILE_UNIT, RESNORM_FILE_UNIT, FILE_NAME_LENGTH, &
       STRING_BUFFER_LENGTH, INTERPOLANT_NAME_LENGTH
  use utils, only: alloc, dealloc, dmsg, DEBUG_LEVEL
  use grid 
  implicit none

  ! process layout
  integer,public :: total_process,total_entries,process_id, &
       left_id,top_id,right_id,bottom_id,back_id,front_id
  character(len=FILE_NAME_LENGTH) :: grid_file_buf
  real, public, dimension(:),allocatable :: left_send_buf,left_recv_buf,&
  top_send_buf,top_recv_buf,right_send_buf,right_recv_buf, &
       bottom_send_buf,bottom_recv_buf,front_send_buf,front_recv_buf,&
       back_send_buf,back_recv_buf
  public :: get_next_token_parallel
  public :: read_layout_file
  public :: get_process_data

contains


  subroutine allocate_buffer_cells()
  implicit none
  integer :: buf
  call dmsg(1, 'parallel', 'allocating buffer cells for MP')  
  !left buffer jmx+1 * kmx + 1 * (dens, x , y, z speeds, pressure)
  
  buf = (jmx+1)*(kmx+1)*5
  call alloc(left_send_buf, 1,buf, &
                    errmsg='Error: Unable to allocate memory for buffer ' // &
                        'variable left_buf.')
  call alloc(right_send_buf, 1,buf, &
                    errmsg='Error: Unable to allocate memory for buffer ' // &
                        'variable right_buf.')
  call alloc(left_recv_buf, 1,buf, &
                    errmsg='Error: Unable to allocate memory for buffer ' // &
                        'variable left_buf.')
  call alloc(right_recv_buf, 1,buf, &
                    errmsg='Error: Unable to allocate memory for buffer ' // &
                        'variable right_buf.')                      
  end subroutine allocate_buffer_cells
  

  subroutine get_process_data()
  implicit none
    ! finds process data
    integer :: ierr
    call MPI_COMM_RANK(MPI_COMM_WORLD,process_id,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,total_process,ierr)

  end subroutine get_process_data

  subroutine get_next_token_parallel(buf)
    !-----------------------------------------------------------
    ! Extract the next token from the layout file
    !
    ! Each token is on a separate line.
    ! There may be multiple comments (lines beginning with #) 
    ! and blank lines in between.
    ! The purpose of this subroutine is to ignore all these 
    ! lines and return the next "useful" line.
    !-----------------------------------------------------------

    implicit none
    character(len=STRING_BUFFER_LENGTH), intent(out) :: buf
    integer :: ios

    do
       read(CONFIG_FILE_UNIT, '(A)', iostat=ios) buf
       if (ios /= 0) then
          print *, 'Error while reading config file.'
          print *, 'Current buffer length is set to: ', &
               STRING_BUFFER_LENGTH
          stop
       end if
       if (index(buf, '#') == 1) then
          ! The current line begins with a hash
          ! Ignore it
          continue
       else if (len_trim(buf) == 0) then
          ! The current line is empty
          ! Ignore it
          continue
       else
          ! A new token has been found
          ! Break out
          exit
       end if
    end do
    call dmsg(0, 'solver', 'get_next_token', 'Returning: ' // trim(buf))

  end subroutine get_next_token_parallel


  subroutine read_layout_file(process_id)
    implicit none
    character(len=FILE_NAME_LENGTH) :: layout_file = "layout.md"
    character(len=STRING_BUFFER_LENGTH) :: buf
    integer,intent(in)::process_id
    integer :: j,i 
    call dmsg(1, 'parallel', 'read_layout_file')

    open(CONFIG_FILE_UNIT, file=layout_file)

    ! Read the parameters from the file
    call get_next_token_parallel(buf)
    read(buf,*)total_process
    call get_next_token_parallel(buf)
    read(buf,*)total_entries
    i = 0
    do while(i < process_id)
       do j = 1, total_entries
          call get_next_token_parallel(buf)
       end do
       i = i+1
    end do
    call get_next_token_parallel(buf) ! process id ignore
    call get_next_token_parallel(buf)
    read(buf,*)left_id
    call get_next_token_parallel(buf)
    read(buf,*)top_id
    call get_next_token_parallel(buf)
    read(buf,*)right_id
    call get_next_token_parallel(buf)
    read(buf,*)bottom_id
    call get_next_token_parallel(buf)
    read(buf,*)back_id
    call get_next_token_parallel(buf)
    read(buf,*)front_id
    call get_next_token_parallel(buf)
    read(buf,*) grid_file_buf         
    !call dmsg(5, 'solver', 'read_config_file', &
    !        msg='scheme_name = ' + scheme_name)
  end subroutine read_layout_file


end module parallel
