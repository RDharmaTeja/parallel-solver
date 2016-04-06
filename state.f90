module state
  !-------------------------------------------------------------------
  ! The state module contains the state variables and the methods that
  ! act on them. 
  !
  ! The state of the system is defined using the density, velocity and
  ! pressure (primitive variables qp) at the grid points. This current
  ! version assumes the grid is in atmost two dimensions. 
  !-------------------------------------------------------------------
  use mpi
  use global, only: FILE_NAME_LENGTH, STATE_FILE_UNIT, OUT_FILE_UNIT, &
       DESCRIPTION_STRING_LENGTH
  use utils, only: alloc, dealloc, dmsg
  use string
  use grid, only: imx, jmx, kmx
  !                   , n_sph_ind, sphere_indices
  use geometry, only: xnx, xny, xnz, ynx, yny, ynz, znx, zny, znz
  use parallel
  implicit none
  private

  ! State variables
  real, public, dimension(:, :, :, :), allocatable, target :: qp
  ! Infinity variables (free stream conditions)
  real, public, dimension(:), allocatable, target :: qp_inf
  ! State variable component aliases
  integer, public :: n_var
  real, public, dimension(:, :, :), pointer :: density
  real, public, dimension(:, :, :), pointer :: x_speed
  real, public, dimension(:, :, :), pointer :: y_speed
  real, public, dimension(:, :, :), pointer :: z_speed
  real, public, dimension(:, :, :), pointer :: pressure
  real, public, pointer :: density_inf
  real, public, pointer :: x_speed_inf
  real, public, pointer :: y_speed_inf
  real, public, pointer :: z_speed_inf
  real, public, pointer :: pressure_inf
  ! Supersonic flag
  logical :: supersonic_flag
  ! Ratio of specific heats (gamma)
  real, public :: gm
  ! Specific gas constant
  real, public :: R_gas
  ! Constants related to viscosity
  real, public :: mu_ref, T_ref, Sutherland_temp, Pr

  ! Public methods
  public :: setup_state
  public :: destroy_state
  public :: set_ghost_cell_data
  !   public :: xi_face_normal_speeds
  !   public :: eta_face_normal_speeds
  !   public :: tau_face_normal_speeds
  public :: writestate, writestate_extra

contains

  subroutine link_aliases()
    implicit none
    call dmsg(1, 'state', 'link_aliases')
    density(0:imx, 0:jmx, 0:kmx) => qp(:, :, :, 1)
    x_speed(0:imx, 0:jmx, 0:kmx) => qp(:, :, :, 2)
    y_speed(0:imx, 0:jmx, 0:kmx) => qp(:, :, :, 3)
    z_speed(0:imx, 0:jmx, 0:kmx) => qp(:, :, :, 4)
    pressure(0:imx, 0:jmx, 0:kmx) => qp(:, :, :, 5)
    density_inf => qp_inf(1)
    x_speed_inf => qp_inf(2)
    y_speed_inf => qp_inf(3)
    z_speed_inf => qp_inf(4)
    pressure_inf => qp_inf(5)
  end subroutine link_aliases

  subroutine unlink_aliases()
    implicit none
    call dmsg(1, 'state', 'unlink_aliases')
    nullify(density)
    nullify(x_speed)
    nullify(y_speed)
    nullify(z_speed)
    nullify(pressure)
    nullify(density_inf)
    nullify(x_speed_inf)
    nullify(y_speed_inf)
    nullify(z_speed_inf)
    nullify(pressure_inf)
  end subroutine unlink_aliases

  subroutine allocate_memory()
    !-----------------------------------------------------------
    ! Allocate memory for the state variables
    !
    ! This assumes that imx and jmx (the grid size) has been set
    ! within the state module.
    !-----------------------------------------------------------

    implicit none

    call dmsg(1, 'state', 'allocate_memory')

    ! The state of the system is defined by the primitive 
    ! variables (density, velocity and pressure) at the grid
    ! cell centers. 
    ! There are (imx - 1) x (jmx - 1) grid cells within the 
    ! domain. We require a row of ghost cells on each boundary.
    ! This current implementation is for a 2D/1D case. 
    call alloc(qp, 0, imx, 0, jmx, 0, kmx, 1, n_var, &
         errmsg='Error: Unable to allocate memory for state ' // &
         'variable qp.')
    call alloc(qp_inf, 1, n_var, &
         errmsg='Error: Unable to allocate memory for state ' // &
         'variable qp_inf.')
  end subroutine allocate_memory

  subroutine deallocate_memory()

    implicit none

    call dmsg(1, 'state', 'deallocate_memory')

    call dealloc(qp)

  end subroutine deallocate_memory

  subroutine setup_state(free_stream_density, free_stream_x_speed, &
       free_stream_y_speed, free_stream_z_speed, &
       free_stream_pressure, state_file)
    !-----------------------------------------------------------
    ! Setup the state module.
    !
    ! This subroutine should be run before the state variables
    ! are initilized. This subroutine allocates the memory for 
    ! state variables and sets up the aliases to refer to the 
    ! components of the state.
    !-----------------------------------------------------------

    implicit none
    real, intent(in) :: free_stream_density
    real, intent(in) :: free_stream_x_speed, free_stream_y_speed, &
         free_stream_z_speed
    real, intent(in) :: free_stream_pressure
    character(len=FILE_NAME_LENGTH), intent(in) :: state_file

    call dmsg(1, 'state', 'setup_state')

    call allocate_memory()
    call link_aliases()
    call init_infinity_values(free_stream_density, &
         free_stream_x_speed, free_stream_y_speed, &
         free_stream_z_speed, free_stream_pressure)
    call set_supersonic_flag()
    call initstate(state_file)

  end subroutine setup_state

  subroutine destroy_state()
    !-----------------------------------------------------------
    ! Destroy the state module.
    !
    ! This subroutine destroys the state module which includes
    ! unlinking the aliases for the state components and 
    ! deallocating the memory held by the state variables. 
    !-----------------------------------------------------------

    implicit none

    call dmsg(1, 'state', 'destroy_state')

    call unlink_aliases()
    call deallocate_memory()

  end subroutine destroy_state
  
  
  subroutine send_recv_interface_lr()
      !-----------------------------------------------------------
    !send and receive data left and right 
    ! i.e along flow direction 
    !-----------------------------------------------------------

    implicit none
    integer :: j,k,count =1
    integer :: ierr,buf
    integer :: status(MPI_STATUS_SIZE)
    if(left_id >= 0 ) then
       !print *, "left id is", process_id,left_id
       ! first send message
   ! left_send_buf(1) = 10;
       count = 1
       do k = 0, kmx
          do j = 0, jmx
             left_send_buf(count) = density(1,j,k)
             count = count+1
          end do
       end do
        do k = 0, kmx
          do j = 0, jmx
             left_send_buf(count) = pressure(1,j,k)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             left_send_buf(count) = x_speed(1,j,k)
             !print *, "left send -  ", process_id ,j,k,count,left_send_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             left_send_buf(count) = y_speed(1,j,k)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             left_send_buf(count) = z_speed(1,j,k)
             count = count+1
          end do
       end do
       !print *,"count is ", count 
       ! send message to left process
       buf = (jmx+1)*(kmx+1)*5       
       !do k = 1, buf
       !print *,'left send - ', process_id, k ,left_send_buf(k)
       !end do       
       call MPI_SEND(left_send_buf,buf,MPI_DOUBLE_PRECISION,left_id,1,MPI_COMM_WORLD, ierr)
       call MPI_RECV(left_recv_buf,buf,MPI_DOUBLE_PRECISION,left_id,1,MPI_COMM_WORLD,status,ierr)
       ! updating solution
       !do k = 1, buf
       !print *,'left recv - ', process_id, k ,left_recv_buf(k)
       !end do       
       count = 1
       do k = 0, kmx
          do j = 0, jmx
              density(0,j,k) = left_recv_buf(count)
               
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             pressure(0,j,k) = left_recv_buf(count)
             
            ! print *, "left recv- ", pressure(0,j,k)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             x_speed(0,j,k) = left_recv_buf(count)
             
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             y_speed(0,j,k) = left_recv_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             z_speed(0,j,k) = left_recv_buf(count)
             count = count+1
          end do
       end do                    
    end if
    
    if(top_id >= 0) then 
       !print *,"top id is ", process_id,top_id
    end if

    if(right_id >= 0) then
       !print *,"right id is ", process_id,right_id
       buf = (jmx+1)*(kmx+1)*5
       call MPI_RECV(right_recv_buf,buf,MPI_DOUBLE_PRECISION,right_id,1,MPI_COMM_WORLD,status,ierr)
       ! updating solution       
       count = 1
       do k = 0, kmx
          do j = 0, jmx
              density(imx,j,k) = right_recv_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             pressure(imx,j,k) = right_recv_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             x_speed(imx,j,k) = right_recv_buf(count)
             !print *, "right recv -  ", process_id ,j,k,count,right_recv_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             y_speed(imx,j,k) = right_recv_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             z_speed(imx,j,k) = right_recv_buf(count)
             count = count+1
          end do
       end do
       ! creating right send  buffer
       count = 1
       do k = 0, kmx
          do j = 0, jmx
             right_send_buf(count) = density(imx-1,j,k)
             
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             right_send_buf(count) = pressure(imx-1,j,k)
             !print *, "left send - ", pressure(imx-1,j,k)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             right_send_buf(count) = x_speed(imx-1,j,k)             
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             right_send_buf(count) = y_speed(imx-1,j,k)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do j = 0, jmx
             right_send_buf(count) = z_speed(imx-1,j,k)
             count = count+1
          end do
       end do
       !do k = 1, buf
       !print *,'right send - ', process_id, k ,right_send_buf(k)
       !end do
       !print *, "total size is", buf 
       call MPI_SEND(right_send_buf,buf,MPI_DOUBLE_PRECISION,right_id,1,MPI_COMM_WORLD, ierr)    
    end if

    if(bottom_id >= 0) then 
       !print *,"bottom id is ", process_id,bottom_id
    end if

    if(back_id >= 0) then
       !print *,"back id is ", process_id,back_id
    end if

    if(front_id >= 0) then
       !print *,"front id is ", process_id,front_id
    end if
    
  end subroutine send_recv_interface_lr


  subroutine send_recv_interface_tb()
      !-----------------------------------------------------------
    !send and receive data top and botom 
    ! i.e along flow direction 
    !-----------------------------------------------------------
    implicit none
    integer :: i,k,count =1
    integer :: ierr,buf
    integer :: status(MPI_STATUS_SIZE)
    call dmsg(1, 'state', 'send recv interface top and bottom')
    if(top_id >= 0 ) then
       !print *, "left id is", process_id,left_id
       ! first send message
   ! left_send_buf(1) = 10;
       count = 1
       do k = 0, kmx
          do i = 0, imx
             top_send_buf(count) = density(i,jmx-1,k)
             count = count+1
          end do
       end do
        do k = 0, kmx
          do i = 0, imx
             top_send_buf(count) = pressure(i,jmx-1,k)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             top_send_buf(count) = x_speed(i,jmx-1,k)
             !print *, "left send -  ", process_id ,j,k,count,left_send_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             top_send_buf(count) = y_speed(i,jmx-1,k)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             top_send_buf(count) = z_speed(i,jmx-1,k)
             count = count+1
          end do
       end do
       !print *,"count is ", count 
       ! send message to left process
       buf = (imx+1)*(kmx+1)*5       
       !do k = 1, buf
       !print *,'left send - ', process_id, k ,left_send_buf(k)
       !end do       
       call MPI_SEND(top_send_buf,buf,MPI_DOUBLE_PRECISION,top_id,1,MPI_COMM_WORLD, ierr)
       call MPI_RECV(top_recv_buf,buf,MPI_DOUBLE_PRECISION,top_id,1,MPI_COMM_WORLD,status,ierr)
       ! updating solution
       !do k = 1, buf
       !print *,'left recv - ', process_id, k ,left_recv_buf(k)
       !end do       
       count = 1
       do k = 0, kmx
          do i = 0, imx
              density(i,jmx,k) = top_recv_buf(count)  
              !print *, "top recv - den ", process_id,top_recv_buf(count)             
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             pressure(i,jmx,k) = top_recv_buf(count)
             !print *, "top recv - pres ",process_id, top_recv_buf(count)
            ! print *, "left recv- ", pressure(0,j,k)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             x_speed(i,jmx,k) = top_recv_buf(count)
             !print *, "top recv - x ", process_id,top_recv_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             y_speed(i,jmx,k) = top_recv_buf(count)
             !print *, "top recv - y ", process_id,top_recv_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             z_speed(i,jmx,k) = top_recv_buf(count)
             !print *, "top recv - x ", process_id, top_recv_buf(count)
             count = count+1
          end do
       end do
       
                           
    end if
    

    if(bottom_id >= 0) then
       !print *,"right id is ", process_id,right_id
       buf = (imx+1)*(kmx+1)*5
       call MPI_RECV(bottom_recv_buf,buf,MPI_DOUBLE_PRECISION,bottom_id,1,MPI_COMM_WORLD,status,ierr)
       ! updating solution 
             
       count = 1
       do k = 0, kmx
          do i = 0, imx
              density(i,0,k) = bottom_recv_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             pressure(i,0,k) = bottom_recv_buf(count)
             
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             x_speed(i,0,k) = bottom_recv_buf(count)
             !print *, "right recv -  ", process_id ,j,k,count,right_recv_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             y_speed(i,0,k) = bottom_recv_buf(count)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             z_speed(i,0,k) = bottom_recv_buf(count)
             count = count+1
          end do
       end do
       ! creating bottom send  buffer
       count = 1
       do k = 0, kmx
          do i = 0, imx
             bottom_send_buf(count) = density(i,1,k)
             
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             bottom_send_buf(count) = pressure(i,1,k)
             !print *, "bottom send - ",process_id, pressure(i,1,k)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             bottom_send_buf(count) = x_speed(i,1,k)             
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             bottom_send_buf(count) = y_speed(i,1,k)
             count = count+1
          end do
       end do
       do k = 0, kmx
          do i = 0, imx
             bottom_send_buf(count) = z_speed(i,1,k)
             count = count+1
          end do
       end do
       !do k = 1, buf
       !print *,'right send - ', process_id, k ,right_send_buf(k)
       !end do
       !print *, "total size is", buf 
       call MPI_SEND(bottom_send_buf,buf,MPI_DOUBLE_PRECISION,bottom_id,1,MPI_COMM_WORLD, ierr)    
    end if
  end subroutine send_recv_interface_tb

  subroutine set_inlet_and_exit_state_variables()
    !-----------------------------------------------------------
    ! Set / extrapolate inlet and exit state variables
    !
    ! The inlet and exit variables are set based on the 
    ! prescribed infinity (free stream) values and extrapolation
    ! of the neighbouring cells. The decision to impose or 
    ! extrapolate is based on the wave speeds and whether or not
    ! the flow is supersonic.
    !-----------------------------------------------------------

    implicit none

    call dmsg(1, 'state', 'set_inlet_and_exit_state_variables')
    
    ! Impose the density, x_speed and y_speed at the inlet
    
    density(0, :, :) = density_inf
    x_speed(0, :, :) = x_speed_inf
    y_speed(0, :, :) = y_speed_inf
    z_speed(0, :, :) = z_speed_inf
    ! Extrapolate these quantities at the exit
    density(imx, :, :) = density(imx - 1, :, :)
    x_speed(imx, :, :) = x_speed(imx - 1, :, :)
    y_speed(imx, :, :) = y_speed(imx - 1, :, :)
    z_speed(imx, :, :) = z_speed(imx - 1, :, :)
    ! If the flow is subsonic, impose the back pressure
    ! Else, impose the inlet pressure
    ! Extrapolate at the other end
    if (supersonic_flag .eqv. .TRUE.) then
       pressure(0, :, :) = pressure_inf
       pressure(imx, :, :) = pressure(imx - 1, :, :)
    else
       pressure(imx, :, :) = pressure_inf
       pressure(0, :, :) = pressure(1, :, :)
    end if
    call send_recv_interface_lr()
    call set_exit_spherical_wall_variables()

  end subroutine set_inlet_and_exit_state_variables

  subroutine set_exit_spherical_wall_variables

    implicit none

    !  integer :: m, i, j, k

    call dmsg(1, 'state', 'set_exit_spherical_wall_variables')

    ! sphere_index contains index of the face at imx. Since the
    ! indices of ghost cells are the same as the face, we use it directly
    ! This subroutine will be called after setting the inlet and exit
    ! So, we are over-writing it
    ! At the same time, flow tangency needs to be applied

    !  do m = 1, n_sph_ind
    !      i = sphere_indices(1, m)
    !      j = sphere_indices(2, m)
    !      k = sphere_indices(3, m)

    !      density(i, j, k) = density(i-1, j, k)
    !      pressure(i, j, k) = pressure(i-1, j, k)
    !      x_speed(i, j, k) = x_speed(i-1, j, k)
    !      y_speed(i, j, k) = y_speed(i-1, j, k)
    !      z_speed(i, j, k) = z_speed(i-1, j, k)

    !      x_speed(i, j, k) = x_speed(i-1, j, k) - (2* &
    !                         ((x_speed(i-1, j, k) * xnx(i, j, k) * xnx(i, j, k)) + &
    !                          (y_speed(i-1, j, k) * xny(i, j, k) * xnx(i, j, k)) + &
    !                          (z_speed(i-1, j, k) * xnz(i, j, k) * xnx(i, j, k))) &
    !                          )
    !      y_speed(i, j, k) = y_speed(i-1, j, k) - (2* &
    !                         ((x_speed(i-1, j, k) * xnx(i, j, k) * xny(i, j, k)) + &
    !                          (y_speed(i-1, j, k) * xny(i, j, k) * xny(i, j, k)) + &
    !                          (z_speed(i-1, j, k) * xnz(i, j, k) * xny(i, j, k))) &
    !                          )
    !      z_speed(i, j, k) = z_speed(i-1, j, k) - (2* &
    !                         ((x_speed(i-1, j, k) * xnx(i, j, k) * xnz(i, j, k)) + &
    !                          (y_speed(i-1, j, k) * xny(i, j, k) * xnz(i, j, k)) + &
    !                          (z_speed(i-1, j, k) * xnz(i, j, k) * xnz(i, j, k))) &
    !                          )
    !  end do    

  end subroutine set_exit_spherical_wall_variables

  subroutine set_front_and_back_ghost_cell_data()
    !-----------------------------------------------------------
    ! Set the state variables for the front and back ghosh cells
    !
    ! The pressure and density for the front and back ghost 
    ! cells is extrapolated from the neighboring interior cells.
    ! The velocity components are computed by applying the 
    ! flow tangency conditions.
    !-----------------------------------------------------------

    implicit none

    call dmsg(1, 'state', 'set_front_and_back_ghost_cell_data')
    
    !if (1 == 2) then
    pressure(:, 0, :) = pressure(:, 1, :)
    density(:, 0, :) = density(:, 1, :)
    pressure(:, jmx, :) = pressure(:, jmx-1, :)    
    density(:, jmx, :) = density(:, jmx-1, :)
    !
    call apply_eta_flow_tangency_conditions()
    !end if

    
  end subroutine set_front_and_back_ghost_cell_data

  subroutine set_top_and_bottom_ghost_cell_data()
    !-----------------------------------------------------------
    ! Set the state variables for the top and bottom ghosh cells
    !
    ! The pressure and density for the top and bottom ghost 
    ! cells is extrapolated from the neighboring interior cells.
    ! The velocity components are computed by applying the 
    ! flow tangency conditions.
    !-----------------------------------------------------------

    implicit none

    call dmsg(1, 'state', 'set_top_and_bottom_ghost_cell_data')

    pressure(:, :,  0) = pressure(:, :, 1)
    pressure(:, :, kmx) = pressure(:, :, kmx-1)
    density(:, :, 0) = density(:, :,  1)
    density(:, :, kmx) = density(:, :, kmx-1)
    call apply_zeta_flow_tangency_conditions()
    
  end subroutine set_top_and_bottom_ghost_cell_data

  subroutine apply_eta_flow_tangency_conditions()
    !-----------------------------------------------------------
    ! Apply the flow tangency conditions for the eta face
    !
    ! The flow tangency conditions ensure that there is no flow
    ! across the boundaries. This is done by ensuring that the
    ! flow is parallel to the boundary.
    !-----------------------------------------------------------

    implicit none

    call dmsg(1, 'state', 'apply_eta_flow_tangency_conditions')

    if (mu_ref .eq. 0.0) then
       ! For the back cells
       x_speed(1:imx-1, jmx, 1:kmx-1) = x_speed(1:imx-1, jmx-1, 1:kmx-1) - &
            (2. * &
            ((x_speed(1:imx-1, jmx-1, 1:kmx-1) * &
            ynx(1:imx-1, jmx, 1:kmx-1) * ynx(1:imx-1, jmx, 1:kmx-1)) &
            + (y_speed(1:imx-1, jmx-1, 1:kmx-1) * &
            yny(1:imx-1, jmx, 1:kmx-1) * ynx(1:imx-1, jmx, 1:kmx-1)) &
            + (z_speed(1:imx-1, jmx-1, 1:kmx-1) * &
            ynz(1:imx-1, jmx, 1:kmx-1) * ynx(1:imx-1, jmx, 1:kmx-1)) &
            ) &
            )
       y_speed(1:imx-1, jmx, 1:kmx-1) = y_speed(1:imx-1, jmx-1, 1:kmx-1) - &
            (2. * &
            ((x_speed(1:imx-1, jmx-1, 1:kmx-1) * &
            ynx(1:imx-1, jmx, 1:kmx-1) * yny(1:imx-1, jmx, 1:kmx-1)) &
            + (y_speed(1:imx-1, jmx-1, 1:kmx-1) * &
            yny(1:imx-1, jmx, 1:kmx-1) * yny(1:imx-1, jmx, 1:kmx-1)) &
            + (z_speed(1:imx-1, jmx-1, 1:kmx-1) * &
            ynz(1:imx-1, jmx, 1:kmx-1) * yny(1:imx-1, jmx, 1:kmx-1)) &
            ) &
            )
       z_speed(1:imx-1, jmx, 1:kmx-1) = z_speed(1:imx-1, jmx-1, 1:kmx-1) - &
            (2. * &
            ((x_speed(1:imx-1, jmx-1, 1:kmx-1) * &
            ynx(1:imx-1, jmx, 1:kmx-1) * ynz(1:imx-1, jmx, 1:kmx-1)) &
            + (y_speed(1:imx-1, jmx-1, 1:kmx-1) * &
            yny(1:imx-1, jmx, 1:kmx-1) * ynz(1:imx-1, jmx, 1:kmx-1)) &
            + (z_speed(1:imx-1, jmx-1, 1:kmx-1) * &
            ynz(1:imx-1, jmx, 1:kmx-1) * ynz(1:imx-1, jmx, 1:kmx-1)) &
            ) &
            )
       ! For the front cells
       x_speed(1:imx-1, 0, 1:kmx-1) = x_speed(1:imx-1, 1, 1:kmx-1) - &
            (2. * &
            ((x_speed(1:imx-1, 1, 1:kmx-1) * &
            ynx(1:imx-1, 1, 1:kmx-1) * ynx(1:imx-1, 1, 1:kmx-1)) &
            + (y_speed(1:imx-1, 1, 1:kmx-1) * &
            yny(1:imx-1, 1, 1:kmx-1) * ynx(1:imx-1, 1, 1:kmx-1)) &
            + (z_speed(1:imx-1, 1, 1:kmx-1) * &
            ynz(1:imx-1, 1, 1:kmx-1) * ynx(1:imx-1, 1, 1:kmx-1)) &
            ) &
            )
       y_speed(1:imx-1, 0, 1:kmx-1) = y_speed(1:imx-1, 1, 1:kmx-1) - &
            (2. * &
            ((x_speed(1:imx-1, 1, 1:kmx-1) * &
            ynx(1:imx-1, 1, 1:kmx-1) * yny(1:imx-1, 1, 1:kmx-1)) &
            + (y_speed(1:imx-1, 1, 1:kmx-1) * &
            yny(1:imx-1, 1, 1:kmx-1) * yny(1:imx-1, 1, 1:kmx-1)) &
            + (z_speed(1:imx-1, 1, 1:kmx-1) * &
            ynz(1:imx-1, 1, 1:kmx-1) * yny(1:imx-1, 1, 1:kmx-1)) &
            ) &
            )
       z_speed(1:imx-1, 0, 1:kmx-1) = z_speed(1:imx-1, 1, 1:kmx-1) - &
            (2. * &
            ((x_speed(1:imx-1, 1, 1:kmx-1) * &
            ynx(1:imx-1, 1, 1:kmx-1) * znx(1:imx-1, 1, 1:kmx-1)) &
            + (y_speed(1:imx-1, 1, 1:kmx-1) * &
            yny(1:imx-1, 1, 1:kmx-1) * znx(1:imx-1, 1, 1:kmx-1)) &
            + (z_speed(1:imx-1, 1, 1:kmx-1) * &
            ynz(1:imx-1, 1, 1:kmx-1) * znx(1:imx-1, 1, 1:kmx-1)) &
            ) &
            )
    else
       x_speed(1:imx-1, jmx, 1:kmx-1) = - x_speed(1:imx-1, jmx-1, 1:kmx-1)
       y_speed(1:imx-1, jmx, 1:kmx-1) = - y_speed(1:imx-1, jmx-1, 1:kmx-1)
       z_speed(1:imx-1, jmx, 1:kmx-1) = - z_speed(1:imx-1, jmx-1, 1:kmx-1)

       x_speed(1:imx-1, 0, 1:kmx-1) = - x_speed(1:imx-1, 1, 1:kmx-1)
       y_speed(1:imx-1, 0, 1:kmx-1) = - y_speed(1:imx-1, 1, 1:kmx-1)
       z_speed(1:imx-1, 0, 1:kmx-1) = - z_speed(1:imx-1, 1, 1:kmx-1)
    end if
    

  end subroutine apply_eta_flow_tangency_conditions

  subroutine apply_zeta_flow_tangency_conditions()
    !-----------------------------------------------------------
    ! Apply the flow tangency conditions
    !
    ! The flow tangency conditions ensure that there is no flow
    ! across the boundaries. This is done by ensuring that the
    ! flow is parallel to the boundary.
    !-----------------------------------------------------------

    implicit none

    call dmsg(1, 'state', 'apply_zeta_flow_tangency_conditions')
    ! For the top cells
    x_speed(1:imx-1, 1:jmx-1, kmx) = x_speed(1:imx-1, 1:jmx-1, kmx-1) - &
         (2. * &
         ((x_speed(1:imx-1, 1:jmx-1, kmx-1) * &
         znx(1:imx-1, 1:jmx-1, kmx) * znx(1:imx-1, 1:jmx-1, kmx)) &
         + (y_speed(1:imx-1, 1:jmx-1, kmx-1) * &
         zny(1:imx-1, 1:jmx-1, kmx) * znx(1:imx-1, 1:jmx-1, kmx)) &
         + (z_speed(1:imx-1, 1:jmx-1, kmx-1) * &
         znz(1:imx-1, 1:jmx-1, kmx) * znx(1:imx-1, 1:jmx-1, kmx)) &
         ) &
         )
    y_speed(1:imx-1, 1:jmx-1, kmx) = y_speed(1:imx-1, 1:jmx-1, kmx-1) - &
         (2. * &
         ((x_speed(1:imx-1, 1:jmx-1, kmx-1) * &
         znx(1:imx-1, 1:jmx-1, kmx) * zny(1:imx-1, 1:jmx-1, kmx)) &
         + (y_speed(1:imx-1, 1:jmx-1, kmx-1) * &
         zny(1:imx-1, 1:jmx-1, kmx) * zny(1:imx-1, 1:jmx-1, kmx)) &
         + (z_speed(1:imx-1, 1:jmx-1, kmx-1) * &
         znz(1:imx-1, 1:jmx-1, kmx) * zny(1:imx-1, 1:jmx-1, kmx)) &
         ) &
         )
    z_speed(1:imx-1, 1:jmx-1, kmx) = z_speed(1:imx-1, 1:jmx-1, kmx-1) - &
         (2. * &
         ((x_speed(1:imx-1, 1:jmx-1, kmx-1) * &
         znx(1:imx-1, 1:jmx-1, kmx) * znz(1:imx-1, 1:jmx-1, kmx)) &
         + (y_speed(1:imx-1, 1:jmx-1, kmx-1) * &
         zny(1:imx-1, 1:jmx-1, kmx) * znz(1:imx-1, 1:jmx-1, kmx)) &
         + (z_speed(1:imx-1, 1:jmx-1, kmx-1) * &
         znz(1:imx-1, 1:jmx-1, kmx) * znz(1:imx-1, 1:jmx-1, kmx)) &
         ) &
         )

    !           x_speed(1:imx-1, 1:jmx-1, kmx) = - x_speed(1:imx-1, 1:jmx-1, kmx-1)
    !           y_speed(1:imx-1, 1:jmx-1, kmx) = - y_speed(1:imx-1, 1:jmx-1, kmx-1)
    !           z_speed(1:imx-1, 1:jmx-1, kmx) = - z_speed(1:imx-1, 1:jmx-1, kmx-1)

    ! For the bottom cells
    x_speed(1:imx-1, 1:jmx-1, 0) = x_speed(1:imx-1, 1:jmx-1, 1) - &
         (2. * &
         ((x_speed(1:imx-1, 1:jmx-1, 1) * &
         znx(1:imx-1, 1:jmx-1, 1) * znx(1:imx-1, 1:jmx-1, 1)) &
         + (y_speed(1:imx-1, 1:jmx-1, 1) * &
         zny(1:imx-1, 1:jmx-1, 1) * znx(1:imx-1, 1:jmx-1, 1)) &
         + (z_speed(1:imx-1, 1:jmx-1, 1) * &
         znz(1:imx-1, 1:jmx-1, 1) * znx(1:imx-1, 1:jmx-1, 1)) &
         ) &
         )
    y_speed(1:imx-1, 1:jmx-1, 0) = y_speed(1:imx-1, 1:jmx-1, 1) - &
         (2. * &
         ((x_speed(1:imx-1, 1:jmx-1, 1) * &
         znx(1:imx-1, 1:jmx-1, 1) * zny(1:imx-1, 1:jmx-1, 1)) &
         + (y_speed(1:imx-1, 1:jmx-1, 1) * &
         zny(1:imx-1, 1:jmx-1, 1) * zny(1:imx-1, 1:jmx-1, 1)) &
         + (z_speed(1:imx-1, 1:jmx-1, 1) * &
         znz(1:imx-1, 1:jmx-1, 1) * zny(1:imx-1, 1:jmx-1, 1)) &
         ) &
         )
    z_speed(1:imx-1, 1:jmx-1, 0) = z_speed(1:imx-1, 1:jmx-1, 1) - &
         (2. * &
         ((x_speed(1:imx-1, 1:jmx-1, 1) * &
         znx(1:imx-1, 1:jmx-1, 1) * znz(1:imx-1, 1:jmx-1, 1)) &
         + (y_speed(1:imx-1, 1:jmx-1, 1) * &
         zny(1:imx-1, 1:jmx-1, 1) * znz(1:imx-1, 1:jmx-1, 1)) &
         + (z_speed(1:imx-1, 1:jmx-1, 1) * &
         znz(1:imx-1, 1:jmx-1, 1) * znz(1:imx-1, 1:jmx-1, 1)) &
         ) &
         )

    !           x_speed(1:imx-1, 1:jmx-1, 0) = - x_speed(1:imx-1, 1:jmx-1, 1)
    !           y_speed(1:imx-1, 1:jmx-1, 0) = - y_speed(1:imx-1, 1:jmx-1, 1)
    !           z_speed(1:imx-1, 1:jmx-1, 0) = - z_speed(1:imx-1, 1:jmx-1, 1)

  end subroutine apply_zeta_flow_tangency_conditions

  subroutine set_ghost_cell_data()
    !-----------------------------------------------------------
    ! Set the data in the ghost cells
    !
    ! The ghost cell data is either imposed or extrapolated from
    ! the neighboring cells inside the domain. The decision to
    ! impose or extrapolate is taken based on certain 
    ! conditions.
    !-----------------------------------------------------------

    implicit none

    call dmsg(1, 'state', 'set_ghost_cell_data')
    call set_inlet_and_exit_state_variables()
    call set_front_and_back_ghost_cell_data()
    call set_top_and_bottom_ghost_cell_data()
    !if (1 == 2) then 
    call send_recv_interface_tb()
    !end if
    
  end subroutine set_ghost_cell_data

  subroutine init_infinity_values(free_stream_density, &
       free_stream_x_speed, free_stream_y_speed, &
       free_stream_z_speed, free_stream_pressure)
    !-----------------------------------------------------------
    ! Set the values of the infinity variables
    !-----------------------------------------------------------

    implicit none
    real, intent(in) :: free_stream_density
    real, intent(in) :: free_stream_x_speed
    real, intent(in) :: free_stream_y_speed
    real, intent(in) :: free_stream_z_speed
    real, intent(in) :: free_stream_pressure

    call dmsg(1, 'state', 'init_infinity_values')

    density_inf = free_stream_density
    x_speed_inf = free_stream_x_speed
    y_speed_inf = free_stream_y_speed
    z_speed_inf = free_stream_z_speed
    pressure_inf = free_stream_pressure

  end subroutine init_infinity_values

  function sound_speed_inf() result(a)
    !-----------------------------------------------------------
    ! Return the free stream speed of sound.
    !-----------------------------------------------------------

    implicit none
    real :: a

    a = sqrt(gm * pressure_inf / density_inf)

  end function sound_speed_inf

  subroutine set_supersonic_flag()
    !-----------------------------------------------------------
    ! Set the supersonic flag based on the infinity conditions.
    !
    ! In the ghost cells, the values of the primitive variables
    ! are set either based on the neighbouring cells' values or 
    ! the infinity values. This splitting is based on the wave
    ! speeds. The supersonic flag is used as an indication as 
    ! to how these cells should get their values. 
    !-----------------------------------------------------------

    implicit none
    real :: avg_inlet_mach

    call dmsg(1, 'state', 'set_supersonic_flag')

    avg_inlet_mach = sqrt(x_speed_inf ** 2. + y_speed_inf ** 2. + &
         z_speed_inf ** 2.) / sound_speed_inf()

    if (avg_inlet_mach >= 1) then
       supersonic_flag = .TRUE.
    else
       supersonic_flag = .FALSE.
    end if

    call dmsg(5, 'state', 'set_supersonic_flag', &
         'Supersonic flag set to ' + supersonic_flag)

  end subroutine set_supersonic_flag

  subroutine initstate(state_file)
    !-----------------------------------------------------------
    ! Initialize the state
    !
    ! If state_file is a tilde (~), then the state should be 
    ! set to the infinity values. Otherwise, read the state_file
    ! to get the state values.
    !-----------------------------------------------------------

    implicit none
    character(len=FILE_NAME_LENGTH), intent(in) :: state_file

    call dmsg(1, 'state', 'initstate')

    if (state_file .eq. '~') then
       ! Set the state to the infinity values
       call init_state_with_infinity_values()
    else
       call readstate(state_file)
    end if

  end subroutine initstate

  subroutine init_state_with_infinity_values()
    !-----------------------------------------------------------
    ! Initialize the state based on the infinity values.
    !-----------------------------------------------------------

    implicit none
    integer :: i

    call dmsg(1, 'state', 'init_state_with_infinity_values')

    do i = 1,n_var
       qp(:, :, :, i) = qp_inf(i)
    end do

  end subroutine init_state_with_infinity_values

  subroutine writestate(outfile, comment)
    !-----------------------------------------------------------
    ! Write the state of the system to a file
    !-----------------------------------------------------------

    implicit none
    character(len=FILE_NAME_LENGTH), intent(in) :: outfile
    character(len=DESCRIPTION_STRING_LENGTH), optional, intent(in) :: &
         comment
    integer :: i, j, k

    call dmsg(1, 'state', 'writestate')

    open(OUT_FILE_UNIT, file=outfile + '.part')

    if (present(comment)) then
       write(OUT_FILE_UNIT, *) & 
            trim('cfd-iitm output (' + comment + ')')
    else
       write(OUT_FILE_UNIT, *) trim('cfd-iitm output')
    end if

    write(OUT_FILE_UNIT, *) 'CELLDATA'
    write(OUT_FILE_UNIT, *) 'Density'
    do k = 1, kmx - 1
       do j = 1, jmx - 1
          do i = 1, imx - 1
             write(OUT_FILE_UNIT, *) density(i, j, k)
          end do
       end do
    end do

    write(OUT_FILE_UNIT, *) 'CELLDATA'
    write(OUT_FILE_UNIT, *) 'Velocity'
    do k = 1, kmx - 1
       do j = 1, jmx - 1
          do i = 1, imx - 1
             write(OUT_FILE_UNIT, *) x_speed(i, j, k), &
                  y_speed(i, j, k), z_speed(i, j, k)
          end do
       end do
    end do

    write(OUT_FILE_UNIT, *) 'CELLDATA'
    write(OUT_FILE_UNIT, *) 'Pressure'
    do k = 1, kmx - 1
       do j = 1, jmx - 1
          do i = 1, imx - 1
             write(OUT_FILE_UNIT, *) pressure(i, j, k)
          end do
       end do
    end do

    close(OUT_FILE_UNIT)

    call rename(outfile + '.part', outfile)

  end subroutine writestate

  subroutine writestate_extra(outfile, comment, extra_var)
    !-----------------------------------------------------------
    ! Write the state of the system to a file
    !-----------------------------------------------------------

    implicit none
    character(len=FILE_NAME_LENGTH), intent(in) :: outfile
    character(len=DESCRIPTION_STRING_LENGTH), optional, intent(in) :: &
         comment
    integer :: i, j, k
    real, dimension(:, :, :) :: extra_var

    call dmsg(5, 'state', 'writestate_extra')

    open(OUT_FILE_UNIT, file=outfile + '.part')

    if (present(comment)) then
       write(OUT_FILE_UNIT, *) & 
            trim('cfd-iitm output (' + comment + ')')
    else
       write(OUT_FILE_UNIT, *) trim('cfd-iitm output')
    end if

    write(OUT_FILE_UNIT, *) 'CELLDATA'
    write(OUT_FILE_UNIT, *) 'Density'
    do k = 1, kmx - 1
       do j = 1, jmx - 1
          do i = 1, imx - 1
             write(OUT_FILE_UNIT, *) density(i, j, k)
          end do
       end do
    end do

    write(OUT_FILE_UNIT, *) 'CELLDATA'
    write(OUT_FILE_UNIT, *) 'Velocity'
    do k = 1, kmx - 1
       do j = 1, jmx - 1
          do i = 1, imx - 1
             write(OUT_FILE_UNIT, *) x_speed(i, j, k), &
                  y_speed(i, j, k), z_speed(i, j, k)
          end do
       end do
    end do

    write(OUT_FILE_UNIT, *) 'CELLDATA'
    write(OUT_FILE_UNIT, *) 'Pressure'
    do k = 1, kmx - 1
       do j = 1, jmx - 1
          do i = 1, imx - 1
             write(OUT_FILE_UNIT, *) pressure(i, j, k)
          end do
       end do
    end do

    write(OUT_FILE_UNIT, *) 'CELLDATA'
    write(OUT_FILE_UNIT, *) 'ExtraVar'
    do k = 1, kmx - 1
       do j = 1, jmx - 1
          do i = 1, imx - 1
             write(OUT_FILE_UNIT, *) extra_var(i, j, k)
          end do
       end do
    end do

    close(OUT_FILE_UNIT)

    call rename(outfile + '.part', outfile)

  end subroutine writestate_extra

  subroutine readstate(state_file)
    !-----------------------------------------------------------
    ! Initialize the state using a state file.
    !
    ! Prior to running this subroutine, the memory for the 
    ! state variables should have been allocated and the 
    ! pointers to alias the components of the state should have 
    ! been associated.
    !-----------------------------------------------------------

    implicit none
    character(len=FILE_NAME_LENGTH), intent(in) :: state_file
    integer :: i, j, k

    call dmsg(1, 'state', 'readstate')

    open(STATE_FILE_UNIT, file=state_file)

    ! Skip the initial comment
    read(STATE_FILE_UNIT, *)

    ! Skip the section header
    read(STATE_FILE_UNIT, *)
    read(STATE_FILE_UNIT, *)
    do k = 1, kmx - 1
       do j = 1, jmx - 1
          do i = 1, imx - 1
             read(STATE_FILE_UNIT, *) density(i, j, k)
          end do
       end do
    end do

    ! Skip the section header
    read(STATE_FILE_UNIT, *)
    read(STATE_FILE_UNIT, *)
    do k = 1, kmx - 1
       do j = 1, jmx - 1
          do i = 1, imx - 1
             read(STATE_FILE_UNIT, *) x_speed(i, j, k), &
                  y_speed(i, j, k), z_speed(i, j, k)
          end do
       end do
    end do

    ! Skip the section header
    read(STATE_FILE_UNIT, *)
    read(STATE_FILE_UNIT, *)
    do k = 1, kmx - 1
       do j = 1, jmx - 1
          do i = 1, imx - 1
             read(STATE_FILE_UNIT, *) pressure(i, j, k)
          end do
       end do
    end do

    close(STATE_FILE_UNIT)

  end subroutine readstate

end module state
