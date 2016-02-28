module solver

    use global, only: CONFIG_FILE_UNIT, RESNORM_FILE_UNIT, FILE_NAME_LENGTH, &
            STRING_BUFFER_LENGTH, INTERPOLANT_NAME_LENGTH
    use utils, only: alloc, dealloc, dmsg, DEBUG_LEVEL
    use string
    use grid, only: imx, jmx, kmx, setup_grid, destroy_grid
!                   sphere_indices, n_sph_ind 
    use geometry, only: xnx, xny, xnz, ynx, yny, ynz, znx, zny, znz, &
            xA, yA, zA, volume, setup_geometry, destroy_geometry
    use state, only: n_var, qp, qp_inf, density, x_speed, y_speed, z_speed, &
            pressure, density_inf, x_speed_inf, y_speed_inf, z_speed_inf, pressure_inf, &
            gm, R_gas, setup_state, destroy_state, set_ghost_cell_data, writestate, &
            mu_ref, T_ref, Sutherland_temp, Pr, writestate_extra
    use face_interpolant, only: interpolant, &
            x_qp_left, x_qp_right, &
            y_qp_left, y_qp_right, &
            z_qp_left, z_qp_right
    use scheme, only: scheme_name, residue, setup_scheme, destroy_scheme, &
            compute_residue

    use viscous, only: gradu_x, gradu_y, gradu_z, gradv_x, gradv_y, gradv_z, &
                       gradw_x, gradw_y, gradw_z, gradT_x, gradT_y, gradT_z
    
    ! parallel module
    use parallel
    implicit none
    private

    real, public :: CFL
    character, public :: time_stepping_method
    real, public :: global_time_step
    character(len=INTERPOLANT_NAME_LENGTH) :: time_step_accuracy
    real, dimension(:, :, :, :), allocatable :: qp_n, dEdx_1, dEdx_2, dEdx_3
    real :: tolerance
    integer, public :: max_iters
    integer, public :: checkpoint_iter
    real, public :: resnorm, resnorm_0
    real, public, dimension(:, :, :), allocatable :: delta_t
    integer, public :: iter
    real :: sim_clock
    real :: speed_inf
    ! buffer results folder
    character(len = 100) :: buffer_folder = 'buffer/output'
    

!   real, dimension(:, :, :, :), allocatable, target :: qp_temp
!   real, dimension(:, :, :), pointer :: density_temp, x_speed_temp, &
!                                        y_speed_temp, z_speed_temp, &
!                                        pressure_temp

    real, dimension(:), allocatable, target :: qp_temp
    real, pointer :: density_temp, x_speed_temp, &
                                         y_speed_temp, z_speed_temp, &
                                         pressure_temp
    real, dimension(:, :, :), pointer :: mass_residue, x_mom_residue, &
                             y_mom_residue, z_mom_residue, energy_residue

    ! Public methods
    public :: setup_solver
    public :: destroy_solver
    public :: step
    public :: converged

    
    contains

        subroutine get_next_token(buf)
            !-----------------------------------------------------------
            ! Extract the next token from the config file
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

        end subroutine get_next_token

        subroutine read_config_file(free_stream_density, &
                free_stream_x_speed, free_stream_y_speed, &
                free_stream_z_speed, free_stream_pressure, &
                grid_file, state_load_file)

            implicit none
            real, intent(out) :: free_stream_density, free_stream_x_speed, &
                    free_stream_y_speed, free_stream_z_speed, free_stream_pressure
            character(len=FILE_NAME_LENGTH), intent(out) :: grid_file
            character(len=FILE_NAME_LENGTH), intent(out) :: state_load_file
            character(len=FILE_NAME_LENGTH) :: config_file = "config.md"
            character(len=STRING_BUFFER_LENGTH) :: buf
            integer :: ios

            call dmsg(1, 'solver', 'read_config_file')
            
            open(CONFIG_FILE_UNIT, file=config_file)

            ! Ignore the config file header
            read(CONFIG_FILE_UNIT, *)
            read(CONFIG_FILE_UNIT, *)
            
            ! Read the parameters from the file

            call get_next_token(buf)
            read(buf, *) scheme_name
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='scheme_name = ' + scheme_name)

            call get_next_token(buf)
            read(buf, *) interpolant
            interpolant = trim(interpolant)
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='interpolant = ' + interpolant)

            call get_next_token(buf)
            read(buf, *) CFL
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='CFL = ' + CFL)

            call get_next_token(buf)
            read(buf, *, iostat=ios) time_stepping_method, global_time_step
            if (ios /= 0) then
                read(buf, *) time_stepping_method
                global_time_step = -1
            end if
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='time_stepping_method = ' + time_stepping_method)
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='global_time_step = ' + global_time_step)

            call get_next_token(buf)
            read(buf, *) time_step_accuracy
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='time_step_accuracy  = ' + time_step_accuracy)

            call get_next_token(buf)
            read(buf, *) tolerance
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='tolerance  = ' + tolerance)

            call get_next_token(buf)
            read(buf, *) grid_file
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='grid_file = ' + grid_file)

            call get_next_token(buf)
            read(buf, *) state_load_file
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='state_load_file = ' + state_load_file)

            call get_next_token(buf)
            read(buf, *) max_iters
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='max_iters = ' + max_iters)

            call get_next_token(buf)
            read(buf, *) checkpoint_iter
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='checkpoint_iter = ' + checkpoint_iter)

            call get_next_token(buf)
            read(buf, *) DEBUG_LEVEL
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='DEBUG_LEVEL = ' + DEBUG_LEVEL)

            call get_next_token(buf)
            read(buf, *) gm
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='gamma = ' + gm)

            call get_next_token(buf)
            read(buf, *) R_gas
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='R_gas = ' + R_gas)
            
            call get_next_token(buf)
            read(buf, *) n_var
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='Number of variables = ' + n_var)

            call get_next_token(buf)
            read(buf, *) free_stream_density
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='free_stream_density = ' + free_stream_density)

            call get_next_token(buf)
            read(buf, *) free_stream_x_speed
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='free_stream_x_speed = ' + free_stream_x_speed)

            call get_next_token(buf)
            read(buf, *) free_stream_y_speed
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='free_stream_y_speed = ' + free_stream_y_speed)

            call get_next_token(buf)
            read(buf, *) free_stream_z_speed
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='free_stream_z_speed = ' + free_stream_z_speed)

            call get_next_token(buf)
            read(buf, *) free_stream_pressure
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='free_stream_pressure = ' + free_stream_pressure)

            call get_next_token(buf)
            read(buf, *) mu_ref
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='mu_reference = ' + mu_ref)

            call get_next_token(buf)
            read(buf, *) T_ref
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='T_reference = ' + T_ref)

            call get_next_token(buf)
            read(buf, *) Sutherland_temp
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='Sutherland temperature = ' + Sutherland_temp)

            call get_next_token(buf)
            read(buf, *) Pr
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='Prandtl Number = ' + Pr)

            close(CONFIG_FILE_UNIT)

        end subroutine read_config_file

        subroutine setup_solver()
            
            implicit none
            real :: free_stream_density
            real :: free_stream_x_speed, free_stream_y_speed, free_stream_z_speed
            real :: free_stream_pressure
            character(len=FILE_NAME_LENGTH) :: grid_file
            character(len=FILE_NAME_LENGTH) :: state_load_file
            
            call dmsg(1, 'solver', 'setup_solver')            
            
            ! parallel calls start
            call get_process_data()
            call read_layout_file(process_id)
            !print *,"process id is",process_id,left_id,top_id,right_id,bottom_id,back_id,front_id,grid_file_name            
			
			! parallel calls ends
            call read_config_file(free_stream_density, free_stream_x_speed, &
                    free_stream_y_speed, free_stream_z_speed, &
                    free_stream_pressure, grid_file, state_load_file)      
            grid_file = grid_file_buf ! changing file name 
            call setup_grid(grid_file)
            call setup_geometry()
            call setup_state(free_stream_density, free_stream_x_speed, &
                    free_stream_y_speed, free_stream_z_speed, &
                    free_stream_pressure, state_load_file)
            call allocate_memory()
            call setup_scheme()
            call link_aliases_solver()
            call initmisc()
            open(RESNORM_FILE_UNIT, file='resnorms')
            call checkpoint()  ! Create an initial dump file

        end subroutine setup_solver

        subroutine destroy_solver()

            implicit none
            
            call dmsg(1, 'solver', 'destroy_solver')

            call destroy_scheme()
            call deallocate_misc()
            call unlink_aliases_solver()
            call destroy_state()
            call destroy_geometry()
            call destroy_grid()
            close(RESNORM_FILE_UNIT)

        end subroutine destroy_solver

        subroutine initmisc()
            
            implicit none
            
            call dmsg(1, 'solver', 'initmisc')

            sim_clock = 0.
            iter = 0
            resnorm = 1.
            resnorm_0 = 1.

        end subroutine initmisc

        subroutine deallocate_misc()

            implicit none
            
            call dmsg(1, 'solver', 'deallocate_misc')

            call dealloc(residue)
            call dealloc(delta_t)

            select case (time_step_accuracy)
                case ("none")
                    ! Do nothing
                    continue
                case ("RK4")
                    call destroy_RK4_time_step()
                case default
                    call dmsg(5, 'solver', 'time_setup_deallocate_memory', &
                                'time step accuracy not recognized.')
                    stop
            end select

        end subroutine deallocate_misc

        subroutine destroy_RK4_time_step()
    
            implicit none

            call dealloc(qp_n)
            call dealloc(dEdx_1)
            call dealloc(dEdx_2)
            call dealloc(dEdx_3)

        end subroutine destroy_RK4_time_step

        subroutine setup_RK4_time_step()
    
            implicit none

            call alloc(qp_n, 0, imx, 0, jmx, 0, kmx, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for qp_n.')
            call alloc(dEdx_1, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for dEdx_1.')
            call alloc(dEdx_2, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for dEdx_2.')
            call alloc(dEdx_3, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for dEdx_3.')

        end subroutine setup_RK4_time_step

        subroutine allocate_memory()

            implicit none
            
            call dmsg(1, 'solver', 'allocate_memory')

            call alloc(delta_t, 1, imx-1, 1, jmx-1, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for delta_t.')
            call alloc(qp_temp, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for qp_temp.')

            select case (time_step_accuracy)
                case ("none")
                    ! Do nothing
                    continue
                case ("RK4")
                    call setup_RK4_time_step()
                case default
                    call dmsg(5, 'solver', 'time_setup_allocate_memory', &
                                'time step accuracy not recognized.')
                    stop
            end select

        end subroutine allocate_memory

        subroutine unlink_aliases_solver()

            implicit none

            nullify(density_temp)
            nullify(x_speed_temp)
            nullify(y_speed_temp)
            nullify(z_speed_temp)
            nullify(pressure_temp)

            nullify(mass_residue)
            nullify(x_mom_residue)
            nullify(y_mom_residue)
            nullify(z_mom_residue)
            nullify(energy_residue)

        end subroutine unlink_aliases_solver

        subroutine link_aliases_solver()

            implicit none

            call dmsg(1, 'solver', 'link_aliases_solver')

        !   density_temp(0:imx, 0:jmx, 0:kmx) => qp_temp(:, :, :, 1)
        !   x_speed_temp(0:imx, 0:jmx, 0:kmx) => qp_temp(:, :, :, 2)
        !   y_speed_temp(0:imx, 0:jmx, 0:kmx) => qp_temp(:, :, :, 3)
        !   z_speed_temp(0:imx, 0:jmx, 0:kmx) => qp_temp(:, :, :, 4)
        !   pressure_temp(0:imx, 0:jmx, 0:kmx) => qp_temp(:, :, :, 5)
            
            density_temp => qp_temp(1)
            x_speed_temp => qp_temp(2)
            y_speed_temp => qp_temp(3)
            z_speed_temp => qp_temp(4)
            pressure_temp => qp_temp(5)
            
            mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 1)
            x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 2)
            y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 3)
            z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 4)
            energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue(:, :, :, 5)

        end subroutine link_aliases_solver

        subroutine compute_local_time_step()
            !-----------------------------------------------------------
            ! Compute the time step to be used at each cell center
            !
            ! Local time stepping can be used to get the solution 
            ! advance towards steady state faster. If only the steady
            ! state solution is required, i.e., transients are 
            ! irrelevant, use local time stepping. 
            !-----------------------------------------------------------

            implicit none

            real :: lmx1, lmx2, lmx3, lmx4, lmx5, lmx6, lmxsum
            real :: x_sound_speed_avg, y_sound_speed_avg, z_sound_speed_avg
            integer :: i, j, k

            call dmsg(1, 'solver', 'compute_local_time_step')

            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
               ! For orientation, refer to the report. The standard i,j,k 
               ! direction are marked. All orientation notations are w.r.t 
               ! to the perspective shown in the image.

               ! Faces with lower index
               x_sound_speed_avg = 0.5 * (sqrt(gm * x_qp_left(i, j, k, 5) / &
                                                    x_qp_left(i, j, k, 1)) + &
                                          sqrt(gm * x_qp_right(i, j, k, 5) / &
                                                    x_qp_right(i, j, k, 1)) )
               y_sound_speed_avg = 0.5 * (sqrt(gm * y_qp_left(i, j, k, 5) / &
                                                    y_qp_left(i, j, k, 1)) + &
                                          sqrt(gm * y_qp_right(i, j, k, 5) / &
                                                    y_qp_right(i, j, k, 1)) )
               z_sound_speed_avg = 0.5 * (sqrt(gm * z_qp_left(i, j, k, 5) / &
                                                    z_qp_left(i, j, k, 1)) + &
                                          sqrt(gm * z_qp_right(i, j, k, 5) / &
                                                    z_qp_right(i, j, k, 1)) )
               
               ! For left face: i.e., lower index face along xi direction
               lmx1 = abs( &
                    (x_speed(i, j, k) * xnx(i, j, k)) + &
                    (y_speed(i, j, k) * xny(i, j, k)) + &
                    (z_speed(i, j, k) * xnz(i, j, k))) + &
                    x_sound_speed_avg
               ! For front face, i.e., lower index face along eta direction
               lmx2 = abs( &
                    (x_speed(i, j, k) * ynx(i, j, k)) + &
                    (y_speed(i, j, k) * yny(i, j, k)) + &
                    (z_speed(i, j, k) * ynz(i, j, k))) + &
                    y_sound_speed_avg
               ! For bottom face, i.e., lower index face along zeta direction
               lmx3 = abs( &
                    (x_speed(i, j, k) * znx(i, j, k)) + &
                    (y_speed(i, j, k) * zny(i, j, k)) + &
                    (z_speed(i, j, k) * znz(i, j, k))) + &
                    z_sound_speed_avg

               ! Faces with higher index
               x_sound_speed_avg = 0.5 * (sqrt(gm * x_qp_left(i+1,j,k,5) / x_qp_left(i+1,j,k,1)) + &
                                          sqrt(gm * x_qp_right(i+1,j,k,5) / x_qp_right(i+1,j,k,1)) )
               y_sound_speed_avg = 0.5 * (sqrt(gm * y_qp_left(i,j+1,k,5) / y_qp_left(i,j+1,k,1)) + &
                                          sqrt(gm * y_qp_right(i,j+1,k,5) / y_qp_right(i,j+1,k,1)) )
               z_sound_speed_avg = 0.5 * (sqrt(gm * z_qp_left(i,j,k+1,5) / z_qp_left(i,j,k+1,1)) + &
                                          sqrt(gm * z_qp_right(i,j,k+1,5) / z_qp_right(i,j,k+1,1)) )
               
               ! For right face, i.e., higher index face along xi direction
               lmx4 = abs( &
                    (x_speed(i+1, j, k) * xnx(i+1, j, k)) + &
                    (y_speed(i+1, j, k) * xny(i+1, j, k)) + &
                    (z_speed(i+1, j, k) * xnz(i+1, j, k))) + &
                    x_sound_speed_avg
               ! For back face, i.e., higher index face along eta direction
               lmx5 = abs( &
                    (x_speed(i, j+1, k) * ynx(i, j+1, k)) + &
                    (y_speed(i, j+1, k) * yny(i, j+1, k)) + &
                    (z_speed(i, j+1, k) * ynz(i, j+1, k))) + &
                    y_sound_speed_avg
               ! For top face, i.e., higher index face along zeta direction
               lmx6 = abs( &
                    (x_speed(i, j, k+1) * znx(i, j, k+1)) + &
                    (y_speed(i, j, k+1) * zny(i, j, k+1)) + &
                    (z_speed(i, j, k+1) * znz(i, j, k+1))) + &
                    z_sound_speed_avg

               lmxsum = (xA(i, j, k) * lmx1) + &
                        (yA(i, j, k) * lmx2) + &
                        (zA(i, j, k) * lmx3) + &
                        (xA(i+1, j, k) * lmx4) + &
                        (yA(i, j+1, k) * lmx5) + &
                        (zA(i, j, k+1) * lmx6)
            
               delta_t(i, j, k) = 1. / lmxsum
               delta_t(i, j, k) = delta_t(i, j, k) * volume(i, j, k) * CFL
              end do
             end do
            end do

        end subroutine compute_local_time_step

        subroutine compute_global_time_step()
            !-----------------------------------------------------------
            ! Compute a common time step to be used at all cell centers
            !
            ! Global time stepping is generally used to get time 
            ! accurate solutions; transients can be studied by 
            ! employing this strategy.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'solver', 'compute_global_time_step')

            if (global_time_step > 0) then
                delta_t = global_time_step
            else
                call compute_local_time_step()
                ! The global time step is the minimum of all the local time
                ! steps.
                delta_t = minval(delta_t)
            end if

        end subroutine compute_global_time_step

        subroutine compute_time_step()
            !-----------------------------------------------------------
            ! Compute the time step to be used
            !
            ! This calls either compute_global_time_step() or 
            ! compute_local_time_step() based on what 
            ! time_stepping_method is set to.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'solver', 'compute_time_step')

            if (time_stepping_method .eq. 'g') then
                call compute_global_time_step()
            else if (time_stepping_method .eq. 'l') then
                call compute_local_time_step()
            else
                call dmsg(5, 'solver', 'compute_time_step', &
                        msg='Value for time_stepping_method (' // &
                            time_stepping_method // ') not recognized.')
                stop
            end if

        end subroutine compute_time_step

        subroutine update_simulation_clock
            !-----------------------------------------------------------
            ! Update the simulation clock
            !
            ! It is sometimes useful to know what the simulation time is
            ! at every iteration so that a comparison with an analytical
            ! solution is possible. Since, the global timesteps used may
            ! not be uniform, we need to track this explicitly.
            !
            ! Of course, it makes sense to track this only if the time 
            ! stepping is global and not local. If the time stepping is
            ! local, the simulation clock is set to -1. If it is global
            ! it is incremented according to the time step found.
            !-----------------------------------------------------------

            implicit none
            if (time_stepping_method .eq. 'g' .and. sim_clock >= 0.) then
                sim_clock = sim_clock + minval(delta_t)
            else if (time_stepping_method .eq. 'l') then
                sim_clock = -1
            end if

        end subroutine update_simulation_clock

        subroutine get_next_solution()

            implicit none

            select case (time_step_accuracy)
                case ("none")
                    call update_solution()
                case ("RK4")
                    call RK4_update_solution()
                case default
                    call dmsg(5, 'solver', 'get_next solution', &
                                'time step accuracy not recognized.')
                    stop
            end select

        end subroutine get_next_solution

        subroutine RK4_update_solution()

            implicit none
            integer :: i, j, k


            ! qp at various stages is not stored but over written
            ! The residue multiplied by the inverse of the jacobian
            ! is stored for the final update equation

            ! Stage 1 is identical to stage (n)
            ! Store qp(n)
            qp_n = qp
            dEdx_1 = get_residue_primitive()
            
            ! Stage 2
            ! Not computing delta_t since qp(1) = qp(n)
            ! Update solution will over write qp
            delta_t = 0.5 * delta_t  ! delta_t(1)
            call update_solution()

            ! Stage 3
            call set_ghost_cell_data()
            call compute_residue()
            call compute_time_step() ! delta_t(2)
            dEdx_2 = get_residue_primitive()
            delta_t = 0.5 * delta_t
            call update_solution()

            ! Stage 4
            call set_ghost_cell_data()
            call compute_residue()
            call compute_time_step() ! delta_t(3)
            dEdx_3 = get_residue_primitive()
            call update_solution()
            ! qp now is qp_4

            ! Use qp(4)
            call set_ghost_cell_data()
            call compute_residue()
            call compute_time_step() ! delta_t(4)

            ! Calculating dEdx_4 in-situ and updating the solution
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                density_temp  = qp_n(i, j, k, 1) - &
                               (((dEdx_1(i, j, k, 1) / 6.0) + &
                                 (dEdx_2(i, j, k, 1) / 3.0) + &
                                 (dEdx_3(i, j, k, 1) / 3.0) + &
                                 (mass_residue(i, j, k) / 6.0)) * &
                                delta_t(i, j, k) / volume(i, j, k))
                x_speed_temp = qp_n(i, j, k, 2) - &
                               (((dEdx_1(i, j, k, 2) / 6.0) + &
                                 (dEdx_2(i, j, k, 2) / 3.0) + &
                                 (dEdx_3(i, j, k, 2) / 3.0) + &
                                 (( (-1 * x_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( x_mom_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
                                ) * delta_t(i, j, k) / volume(i, j, k))
                y_speed_temp = qp_n(i, j, k, 3) - &
                               (((dEdx_1(i, j, k, 3) / 6.0) + &
                                 (dEdx_2(i, j, k, 3) / 3.0) + &
                                 (dEdx_3(i, j, k, 3) / 3.0) + &
                                 (( (-1 * y_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( y_mom_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
                                ) * delta_t(i, j, k) / volume(i, j, k))
                z_speed_temp = qp_n(i, j, k, 4) - &
                               (((dEdx_1(i, j, k, 4) / 6.0) + &
                                 (dEdx_2(i, j, k, 4) / 3.0) + &
                                 (dEdx_3(i, j, k, 4) / 3.0) + &
                                 (( (-1 * z_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( z_mom_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
                                ) * delta_t(i, j, k) / volume(i, j, k))
                pressure_temp = qp_n(i, j, k, 5) - &
                               (((dEdx_1(i, j, k, 5) / 6.0) + &
                                 (dEdx_2(i, j, k, 5) / 3.0) + &
                                 (dEdx_3(i, j, k, 5) / 3.0) + &
                                 (( (0.5 * (gm - 1.) * ( x_speed(i, j, k) ** 2. + &
                                                         y_speed(i, j, k) ** 2. + &
                                                         z_speed(i, j, k) ** 2.) * &
                                                        mass_residue(i, j, k)) + &
                       (- (gm - 1.) * x_speed(i, j, k) * x_mom_residue(i, j, k)) + &
                       (- (gm - 1.) * y_speed(i, j, k) * y_mom_residue(i, j, k)) + &
                       (- (gm - 1.) * z_speed(i, j, k) * z_mom_residue(i, j, k)) + &
                       ((gm - 1.) * energy_residue(i, j, k)) ) / 6.0) &
                                ) * delta_t(i, j, k) / volume(i, j, k))
            
                density(i, j, k) = density_temp
                x_speed(i, j, k) = x_speed_temp
                y_speed(i, j, k) = y_speed_temp
                z_speed(i, j, k) = z_speed_temp
                pressure(i, j, k) = pressure_temp
              end do
             end do
            end do

            if (any(density < 0) .or. any(pressure < 0)) then
                call dmsg(5, 'solver', 'update_solution', &
                        'ERROR: Some density or pressure is negative.')
            end if

        end subroutine RK4_update_solution

        function get_residue_primitive() result(dEdx)

            implicit none

            real, dimension(1:imx-1, 1:jmx-1, 1:kmx-1, n_var) :: dEdx
      !     integer, dimension(1:4) :: dEsize
      !     integer, dimension(1:3) :: ressize

      !     dEsize = shape(dEdx)
      !     ressize = shape(mass_residue)

      !     print *, dEsize(1), dEsize(2), dEsize(3), dEsize(4)
      !     print *, ressize(1), ressize(2), ressize(3)
      !     (1:imx-1, 1:jmx-1, 1:kmx-1)
            dEdx(:, :, :, 1) = mass_residue
            dEdx(:, :, :, 2) = ( (-1 * x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                                     mass_residue) + &
                             ( x_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
            dEdx(:, :, :, 3) = ( (-1 * y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                                     mass_residue) + &
                             ( y_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
            dEdx(:, :, :, 4) = ( (-1 * z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                                     mass_residue) + &
                             ( z_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
            dEdx(:, :, :, 5) = ( (0.5 * (gm - 1.) * ( x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2. + &
                                                      y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2. + &
                                                      z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2.) * &
                                                    mass_residue) + &
                       (- (gm - 1.) * x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * x_mom_residue) + &
                       (- (gm - 1.) * y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * y_mom_residue) + &
                       (- (gm - 1.) * z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * z_mom_residue) + &
                       ((gm - 1.) * energy_residue) )

        end function get_residue_primitive

        subroutine update_solution()
            !-----------------------------------------------------------
            ! Update the solution using the residue and time step
            !-----------------------------------------------------------

            implicit none
            integer :: i, j, k
            
            call dmsg(1, 'solver', 'update_solution')

            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
               density_temp = density(i, j, k) - &
                            (mass_residue(i, j, k) * &
                            delta_t(i, j, k) / volume(i, j, k))

               x_speed_temp = x_speed(i, j, k) - &
                            (( (-1 * x_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( x_mom_residue(i, j, k) / density(i, j, k)) ) * &
                            delta_t(i, j, k) / volume(i, j, k))

               y_speed_temp = y_speed(i, j, k) - &
                            (( (-1 * y_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( y_mom_residue(i, j, k) / density(i, j, k)) ) * &
                            delta_t(i, j, k) / volume(i, j, k))

               z_speed_temp = z_speed(i, j, k) - &
                            (( (-1 * z_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( z_mom_residue(i, j, k) / density(i, j, k)) ) * &
                            delta_t(i, j, k) / volume(i, j, k))

               pressure_temp = pressure(i, j, k) - &
                   ( ( (0.5 * (gm - 1.) * ( x_speed(i, j, k) ** 2. + &
                                            y_speed(i, j, k) ** 2. + &
                                            z_speed(i, j, k) ** 2.) * &
                                          mass_residue(i, j, k)) + &
                       (- (gm - 1.) * x_speed(i, j, k) * x_mom_residue(i, j, k)) + &
                       (- (gm - 1.) * y_speed(i, j, k) * y_mom_residue(i, j, k)) + &
                       (- (gm - 1.) * z_speed(i, j, k) * z_mom_residue(i, j, k)) + &
                       ((gm - 1.) * energy_residue(i, j, k)) ) * &
                       delta_t(i, j, k) / volume(i, j, k) ) 

               density(i, j, k) = density_temp
               x_speed(i, j, k) = x_speed_temp
               y_speed(i, j, k) = y_speed_temp
               z_speed(i, j, k) = z_speed_temp
               pressure(i, j, k) = pressure_temp
              end do
             end do
            end do
            
       !    density(1:imx-1, 1:jmx-1, 1:kmx-1) = density_temp(1:imx-1, 1:jmx-1, 1:kmx-1)
       !    x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) = x_speed_temp(1:imx-1, 1:jmx-1, 1:kmx-1)
       !    y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) = y_speed_temp(1:imx-1, 1:jmx-1, 1:kmx-1)
       !    z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) = z_speed_temp(1:imx-1, 1:jmx-1, 1:kmx-1)
       !    pressure(1:imx-1, 1:jmx-1, 1:kmx-1) = pressure_temp(1:imx-1, 1:jmx-1, 1:kmx-1)

         

         !  qp_temp(1:imx-1, 1:jmx-1, 1:kmx-1,  1) = qp(1:imx-1, 1:jmx-1, 1:kmx-1, 1) - &
         !          (residue(:, :, :, 1) * &
         !          delta_t(:, :, :) / volume(:, :, :))
         !  qp_temp(1:imx-1, 1:jmx-1, 1:kmx-1, 2) = qp(1:imx-1, 1:jmx-1, 1:kmx-1, 2) - &
         !          (( (-1. * qp(1:imx-1, 1:jmx-1, 1:kmx-1, 2) / qp(1:imx-1, 1:jmx-1, 1:kmx-1, 1) * residue(:, :, :, 1)) + &
         !             ( residue(:, :, :, 2) / qp(1:imx-1, 1:jmx-1, 1:kmx-1, 1) )) * &
         !          delta_t(:, :, :) / volume(:, :, :))
         !  qp_temp(1:imx-1, 1:jmx-1, 1:kmx-1, 3) = qp(1:imx-1, 1:jmx-1, 1:kmx-1, 3) - &
         !          (( (-1. * qp(1:imx-1, 1:jmx-1, 1:kmx-1, 3) / qp(1:imx-1, 1:jmx-1, 1:kmx-1, 1) * residue(:, :, :, 1)) + &
         !             ( residue(:, :, :, 3) / qp(1:imx-1, 1:jmx-1, 1:kmx-1, 1) )) * &
         !          delta_t(:, :, :) / volume(:, :, :))
         !  qp_temp(1:imx-1, 1:jmx-1, 1:kmx-1, 4) = qp(1:imx-1, 1:jmx-1, 1:kmx-1, 4) - &
         !          (( (-1. * qp(1:imx-1, 1:jmx-1, 1:kmx-1, 4) / qp(1:imx-1, 1:jmx-1, 1:kmx-1, 1) * residue(:, :, :, 1)) + &
         !             ( residue(:, :, :, 4) / qp(1:imx-1, 1:jmx-1, 1:kmx-1, 1) )) * &
         !          delta_t(:, :, :) / volume(:, :, :))
         !  qp_temp(1:imx-1, 1:jmx-1, 1:kmx-1, 5) = qp(1:imx-1, 1:jmx-1, 1:kmx-1, 5) - &
         !          ( ( (0.5 * (gm - 1.) * ( qp(1:imx-1, 1:jmx-1, 1:kmx-1, 2)**2. + qp(1:imx-1, 1:jmx-1, 1:kmx-1, 3)**2. + &
         !              qp(1:imx-1, 1:jmx-1, 1:kmx-1, 4)**2.) * residue(:, :, :, 1)) + &
         !             (- (gm - 1.) * qp(1:imx-1, 1:jmx-1, 1:kmx-1, 2) * residue(:, :, :, 2)) + &
         !             (- (gm - 1.) * qp(1:imx-1, 1:jmx-1, 1:kmx-1, 3) * residue(:, :, :, 3)) + &
         !             (- (gm - 1.) * qp(1:imx-1, 1:jmx-1, 1:kmx-1, 4) * residue(:, :, :, 4)) + &
         !             ((gm - 1.) * residue(:, :, :, 5)) ) * &
         !          delta_t(:, :, :) / volume(:, :, :) )


         !  qp(1:imx-1, 1:jmx-1, 1:kmx-1, :) = qp_temp(1:imx-1, 1:jmx-1, 1:kmx-1, :)
         
            if (any(density < 0) .or. any(pressure < 0)) then
                call dmsg(5, 'solver', 'update_solution', &
                        'ERROR: Some density or pressure is negative.')
                stop
            end if

        end subroutine update_solution

        subroutine checkpoint()
            !-----------------------------------------------------------
            ! Create a checkpoint dump file if the time has come
            !-----------------------------------------------------------

            implicit none

            character(len=FILE_NAME_LENGTH) :: filename

            if (checkpoint_iter .ne. 0) then
                if (mod(iter, checkpoint_iter) == 0) then
                    write(filename, '(A,I2.2,A,I5.5,A)') 'buffer/process_',process_id,'/output', iter, '.fvtk'
                    print *, filename
                    call writestate(filename, 'Simulation clock: ' + sim_clock)
                !   call writestate_extra(filename, 'Simulation clock: ' + sim_clock, gradu_y(1:imx-1, 1:jmx-1, 1:kmx-1))
                    call dmsg(3, 'solver', 'checkpoint', &
                            'Checkpoint created at iteration: ' + iter)
                end if
            end if

        end subroutine checkpoint

        subroutine step()
            !-----------------------------------------------------------
            ! Perform one time step iteration
            !
            ! This subroutine performs one iteration by stepping through
            ! time once.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'solver', 'step')

            call set_ghost_cell_data()
            call compute_residue()
            call dmsg(1, 'solver', 'step', 'Residue computed.')
            call compute_time_step()
         !  call update_solution()
            call get_next_solution()
            call update_simulation_clock()
            iter = iter + 1
            call compute_residue_norm()
            if (iter .eq. 1) then
                resnorm_0 = resnorm
            end if
            write(RESNORM_FILE_UNIT, *) resnorm
            call checkpoint()

        end subroutine step

        subroutine compute_residue_norm()

            implicit none
            
            call dmsg(1, 'solver', 'compute_residue_norm')

            speed_inf = sqrt(x_speed_inf ** 2. + y_speed_inf ** 2. + &
                             z_speed_inf ** 2.)
            
            resnorm = sum(sqrt( &
                    (residue(:, :, :, 1) / &
                        (density_inf * speed_inf)) ** 2. + &
                    (residue(:, :, :, 2) / &
                        (density_inf * speed_inf ** 2.)) ** 2. + &
                    (residue(:, :, :, 3) / &
                        (density_inf * speed_inf ** 2.)) ** 2. + &
                    (residue(:, :, :, 4) / &
                        (density_inf * speed_inf ** 2.)) ** 2. + &
                    (residue(:, :, :, 5) / &
                        (density_inf * speed_inf ** 3.)) ** 2. &
                    ))

        end subroutine compute_residue_norm

        function converged() result(c)
            !-----------------------------------------------------------
            ! Check if the solution seems to have converged
            !
            ! The solution is said to have converged if the change in 
            ! the residue norm is "negligible".
            !-----------------------------------------------------------

            implicit none
            logical :: c
            
            call dmsg(1, 'solver', 'converged')

            if (resnorm / resnorm_0 < tolerance) then
                c = .TRUE.
            end if
            c = .FALSE.

        end function converged

end module solver

