module van_leer
    !-------------------------------------------------------------------
    ! The Van-Leer scheme is a type of flux-splitting scheme
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx, kmx
!                   sphere_indices, n_sph_ind
    use geometry, only: xnx, xny, xnz, ynx, yny, ynz, znx, zny, znz, xA, yA, zA
    use state, only: gm, n_var
    use face_interpolant, only: x_qp_left, x_qp_right, y_qp_left, y_qp_right, &
                z_qp_left, z_qp_right, &
            x_density_left, x_x_speed_left, x_y_speed_left, x_z_speed_left, &
                x_pressure_left, &
            x_density_right, x_x_speed_right, x_y_speed_right, x_z_speed_right, &
                x_pressure_right, &
            y_density_left, y_x_speed_left, y_y_speed_left, y_z_speed_left, &
                y_pressure_left, &
            y_density_right, y_x_speed_right, y_y_speed_right, y_z_speed_right, &
                y_pressure_right, &
            z_density_left, z_x_speed_left, z_y_speed_left, z_z_speed_left, &
                z_pressure_left, &
            z_density_right, z_x_speed_right, z_y_speed_right, z_z_speed_right, &
                z_pressure_right
 !  use scheme, only: F, G, H, residue

    implicit none
    private

    real, public, dimension(:, :, :, :), allocatable, target :: F, G, H

    !TODO: Integrate here with boundary condition module
    
   !real, dimension(:, :, :), allocatable :: x_sound_speed_avg
   !real, dimension(:, :, :), allocatable :: x_M_perp_left, x_M_perp_right
   !real, dimension(:, :, :), allocatable :: y_sound_speed_avg
   !real, dimension(:, :, :), allocatable :: y_M_perp_left, y_M_perp_right
   !real, dimension(:, :, :), allocatable :: z_sound_speed_avg
   !real, dimension(:, :, :), allocatable :: z_M_perp_left, z_M_perp_right
   !real, dimension(:, :, :), allocatable :: x_alpha_plus, x_alpha_minus
   !real, dimension(:, :, :), allocatable :: y_alpha_plus, y_alpha_minus
   !real, dimension(:, :, :), allocatable :: z_alpha_plus, z_alpha_minus
   !real, dimension(:, :, :), allocatable :: x_beta_left, x_beta_right
   !real, dimension(:, :, :), allocatable :: y_beta_left, y_beta_right
   !real, dimension(:, :, :), allocatable :: z_beta_left, z_beta_right
   !real, dimension(:, :, :), allocatable :: x_c_plus, x_c_minus
   !real, dimension(:, :, :), allocatable :: y_c_plus, y_c_minus
   !real, dimension(:, :, :), allocatable :: z_c_plus, z_c_minus
   !real, dimension(:, :, :), allocatable :: x_scrD_plus, x_scrD_minus
   !real, dimension(:, :, :), allocatable :: y_scrD_plus, y_scrD_minus
   !real, dimension(:, :, :), allocatable :: z_scrD_plus, z_scrD_minus

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_fluxes
    public :: get_residue

    ! Public members which can be used in schemes derived from Van Leer
!   public :: compute_residue

!   public :: compute_xi_face_quantities
!   public :: compute_eta_face_quantities
!   public :: compute_tau_face_quantities
!   public :: x_M_perp_left, x_M_perp_right
!   public :: y_M_perp_left, y_M_perp_right
!   public :: z_M_perp_left, z_M_perp_right
!   public :: x_beta_left, x_beta_right
!   public :: y_beta_left, y_beta_right
!   public :: z_beta_left, z_beta_right
!   public :: x_c_plus, x_c_minus
!   public :: y_c_plus, y_c_minus
!   public :: z_c_plus, z_c_minus

    contains

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'van_leer', 'setup_scheme')

            call alloc(F, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'F - van_leer.')
            call alloc(G, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'G - van_leer.')
            call alloc(H, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'H - van_leer.')

!           call alloc(x_M_perp_left, 1, imx, 1, jmx-1, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'x_M_perp_left.')
!           call alloc(x_M_perp_right, 1, imx, 1, jmx-1, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'x_M_perp_right.')
!           call alloc(y_M_perp_left, 1, imx-1, 1, jmx, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'y_M_perp_left.')
!           call alloc(y_M_perp_right, 1, imx-1, 1, jmx, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'y_M_perp_right.')
!           call alloc(z_M_perp_left, 1, imx-1, 1, jmx-1 , 1, kmx, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'z_M_perp_left.')
!           call alloc(z_M_perp_right, 1, imx-1, 1, jmx-1, 1, kmx, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'z_M_perp_right.')

!           call alloc(x_alpha_plus, 1, imx, 1, jmx-1, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'x_alpha_plus.')
!           call alloc(x_alpha_minus, 1, imx, 1, jmx-1, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'x_alpha_minus.')
!           call alloc(y_alpha_plus, 1, imx-1, 1, jmx, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'y_alpha_plus.')
!           call alloc(y_alpha_minus, 1, imx-1, 1, jmx, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'y_alpha_minus.')
!           call alloc(z_alpha_plus, 1, imx-1, 1, jmx-1, 1, kmx, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'z_alpha_plus.')
!           call alloc(z_alpha_minus, 1, imx-1, 1, jmx-1, 1, kmx, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'z_alpha_minus.')

!           call alloc(x_beta_left, 1, imx, 1, jmx-1, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for x_beta_left.')
!           call alloc(x_beta_right, 1, imx, 1, jmx-1, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'x_beta_right.')
!           call alloc(y_beta_left, 1, imx-1, 1, jmx, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for y_beta_left.')
!           call alloc(y_beta_right, 1, imx-1, 1, jmx, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'y_beta_right.')
!           call alloc(z_beta_left, 1, imx-1, 1, jmx-1, 1, kmx, &
!                   errmsg='Error: Unable to allocate memory for z_beta_left.')
!           call alloc(z_beta_right, 1, imx-1, 1, jmx-1, 1, kmx, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'z_beta_right.')

!           call alloc(x_c_plus, 1, imx, 1, jmx-1, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for x_c_plus.')
!           call alloc(x_c_minus, 1, imx, 1, jmx-1, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for x_c_minus.')
!           call alloc(y_c_plus, 1, imx-1, 1, jmx, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for y_c_plus.')
!           call alloc(y_c_minus, 1, imx-1, 1, jmx, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for y_c_minus.')
!           call alloc(z_c_plus, 1, imx-1, 1, jmx-1, 1, kmx, &
!                   errmsg='Error: Unable to allocate memory for z_c_plus.')
!           call alloc(z_c_minus, 1, imx-1, 1, jmx-1, 1, kmx, &
!                   errmsg='Error: Unable to allocate memory for z_c_minus.')

!           call alloc(x_scrD_plus, 1, imx, 1, jmx-1, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for x_scrD_plus.')
!           call alloc(x_scrD_minus, 1, imx, 1, jmx-1, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'x_scrD_minus.')
!           call alloc(y_scrD_plus, 1, imx-1, 1, jmx, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for y_scrD_plus.')
!           call alloc(y_scrD_minus, 1, imx-1, 1, jmx, 1, kmx-1, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'y_scrD_minus.')
!           call alloc(z_scrD_plus, 1, imx-1, 1, jmx-1, 1, kmx, &
!                   errmsg='Error: Unable to allocate memory for z_scrD_plus.')
!           call alloc(z_scrD_minus, 1, imx-1, 1, jmx-1, 1, kmx, &
!                   errmsg='Error: Unable to allocate memory for ' // &
!                       'z_scrD_minus.')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'van_leer', 'destroy_scheme')
            
            call dealloc(F)
            call dealloc(G)
            call dealloc(H)

!           call dealloc(x_M_perp_left)
!           call dealloc(x_M_perp_right)
!           call dealloc(y_M_perp_left)
!           call dealloc(y_M_perp_right)
!           call dealloc(z_M_perp_left)
!           call dealloc(z_M_perp_right)

!           call dealloc(x_alpha_plus)
!           call dealloc(x_alpha_minus)
!           call dealloc(y_alpha_plus)
!           call dealloc(y_alpha_minus)
!           call dealloc(z_alpha_plus)
!           call dealloc(z_alpha_minus)

!           call dealloc(x_beta_left)
!           call dealloc(x_beta_right)
!           call dealloc(y_beta_left)
!           call dealloc(y_beta_right)
!           call dealloc(z_beta_left)
!           call dealloc(z_beta_right)

!           call dealloc(x_c_plus)
!           call dealloc(x_c_minus)
!           call dealloc(y_c_plus)
!           call dealloc(y_c_minus)
!           call dealloc(z_c_plus)
!           call dealloc(z_c_minus)

!           call dealloc(x_scrD_plus)
!           call dealloc(x_scrD_minus)
!           call dealloc(y_scrD_plus)
!           call dealloc(y_scrD_minus)
!           call dealloc(z_scrD_plus)
!           call dealloc(z_scrD_minus)

        end subroutine destroy_scheme

!       function check_if_sphere_wall(i, j, k) result(res)
!           !-----------------------------------------------------------
!           ! Checks if given indices are part of the spherical wall 
!           !-----------------------------------------------------------
!           integer, intent(in) :: i, j, k
!           logical res

!           integer :: m

!           res = .FALSE.

!           do m = 1, n_sph_ind
!               if ((sphere_indices(1, m) .eq. i) .and. (sphere_indices(2, m) .eq. j) &
!                   .and.(sphere_indices(3, m) .eq. k)) then
!                   res = .TRUE.
!                   return
!               end if
!           end do
!       
!       end function check_if_sphere_wall

        subroutine compute_F_flux()
            !-----------------------------------------------------------
            ! Compute the flux 'F' along the xi faces
            !-----------------------------------------------------------
            
            implicit none

            real, dimension(1:n_var) :: F_plus, F_minus
            real :: M_perp_left, M_perp_right
            real :: alpha_plus, alpha_minus
            real :: beta_left, beta_right
            real :: M_plus, M_minus
            real :: D_plus, D_minus
            real :: c_plus, c_minus
            real :: scrD_plus, scrD_minus
            real :: sound_speed_avg, face_normal_speeds
           !real :: sound_speed_left, sound_speed_right
            integer :: i , j, k
            
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx
                sound_speed_avg = 0.5 * (sqrt(gm * x_pressure_left(i, j, k) / &
                                            x_density_left(i, j, k) ) + &
                                          sqrt(gm * x_pressure_right(i, j, k) / &
                                            x_density_right(i, j, k) ) )
                
                ! Compute '+' direction quantities

                face_normal_speeds = x_x_speed_left(i, j, k) * xnx(i, j, k) + &
                                     x_y_speed_left(i, j, k) * xny(i, j, k) + &
                                     x_z_speed_left(i, j, k) * xnz(i, j, k)
                M_perp_left = face_normal_speeds / sound_speed_avg
                alpha_plus = 0.5 * (1.0 + sign(1.0, M_perp_left))
                beta_left = -max(0, 1 - floor(abs(M_perp_left)))
                M_plus = 0.25 * ((1. + M_perp_left) ** 2.)
                D_plus = 0.25 * ((1. + M_perp_left) ** 2.) * (2. - M_perp_left)
                c_plus = (alpha_plus * (1.0 + beta_left) * M_perp_left) - &
                          beta_left * M_plus
                scrD_plus = (alpha_plus * (1. + beta_left)) - &
                        (beta_left * D_plus)

                ! First construct the F mass flux
                F_plus(1) = x_density_left(i, j, k) * sound_speed_avg * c_plus
                
                ! Is convective flux zero anywhere?
!               if (check_if_sphere_wall(i, j, k)) then
!                   F_plus(1) = 0.
!               end if


                ! Construct other fluxes in terms of the F mass flux
                F_plus(2) = (F_plus(1) * x_x_speed_left(i, j, k)) + &
                            (scrD_plus * x_pressure_left(i, j, k) * xnx(i, j, k))
                F_plus(3) = (F_plus(1) * x_y_speed_left(i, j, k)) + &
                            (scrD_plus * x_pressure_left(i, j, k) * xny(i, j, k))
                F_plus(4) = (F_plus(1) * x_z_speed_left(i, j, k)) + &
                            (scrD_plus * x_pressure_left(i, j, k) * xnz(i, j, k))
                F_plus(5) = F_plus(1) * &
                            ((0.5 * (x_x_speed_left(i, j, k) ** 2. + &
                                     x_y_speed_left(i, j, k) ** 2. + &
                                     x_z_speed_left(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * x_pressure_left(i, j, k) / &
                             x_density_left(i, j, k)))

                ! Multiply in the face areas
                F_plus(1) = F_plus(1) * xA(i, j, k)
                F_plus(2) = F_plus(2) * xA(i, j, k)
                F_plus(3) = F_plus(3) * xA(i, j, k)
                F_plus(4) = F_plus(4) * xA(i, j, k)
                F_plus(5) = F_plus(5) * xA(i, j, k)

                ! Compute '-' direction quantities

                face_normal_speeds = x_x_speed_right(i, j, k) * xnx(i, j, k) + &
                                     x_y_speed_right(i, j, k) * xny(i, j, k) + &
                                     x_z_speed_right(i, j, k) * xnz(i, j, k)
                M_perp_right = face_normal_speeds / sound_speed_avg
                alpha_minus = 0.5 * (1.0 - sign(1.0, M_perp_right))
                beta_right = -max(0, 1 - floor(abs(M_perp_right)))
                M_minus = -0.25 * ((1. - M_perp_right) ** 2.)
                D_minus = 0.25 * ((1. - M_perp_right) ** 2.) * (2. + M_perp_right)
                c_minus = (alpha_minus * (1.0 + beta_right) * M_perp_right) - &
                          beta_right * M_minus
                scrD_minus = (alpha_minus * (1. + beta_right)) - &
                             (beta_right * D_minus)

                ! First construct the F mass flux
                F_minus(1) = x_density_right(i, j, k) * sound_speed_avg * c_minus
                
!               ! Is convective flux zero anywhere?
!               if (check_if_sphere_wall(i, j, k)) then
!                   F_minus(1) = 0.
!               end if

                ! Construct other fluxes in terms of the F mass flux
                F_minus(2) = (F_minus(1) * x_x_speed_right(i, j, k)) + &
                             (scrD_minus * x_pressure_right(i, j, k) * xnx(i, j, k))
                F_minus(3) = (F_minus(1) * x_y_speed_right(i, j, k)) + &
                             (scrD_minus * x_pressure_right(i, j, k) * xny(i, j, k))
                F_minus(4) = (F_minus(1) * x_z_speed_right(i, j, k)) + &
                             (scrD_minus * x_pressure_right(i, j, k) * xnz(i, j, k))
                F_minus(5) = F_minus(1) * &
                            ((0.5 * (x_x_speed_right(i, j, k) ** 2. + &
                                     x_y_speed_right(i, j, k) ** 2. + &
                                     x_z_speed_right(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * x_pressure_right(i, j, k) / &
                             x_density_right(i, j, k)))
         
                ! Multiply in the face areas
                F_minus(1) = F_minus(1) * xA(i, j, k)
                F_minus(2) = F_minus(2) * xA(i, j, k)
                F_minus(3) = F_minus(3) * xA(i, j, k)
                F_minus(4) = F_minus(4) * xA(i, j, k)
                F_minus(5) = F_minus(5) * xA(i, j, k)

                ! Get the total flux for a face
                F(i, j, k, :) = F_plus(:) + F_minus(:)
              end do
             end do
            end do 

        end subroutine compute_F_flux

        subroutine compute_G_flux()
            !-----------------------------------------------------------
            ! Compute the flux 'G' along the eta faces
            !-----------------------------------------------------------
            
            implicit none

            real, dimension(1:n_var) :: G_plus, G_minus
            real :: M_perp_left, M_perp_right
            real :: alpha_plus, alpha_minus
            real :: beta_left, beta_right
            real :: M_plus, M_minus
            real :: D_plus, D_minus
            real :: c_plus, c_minus
            real :: scrD_plus, scrD_minus
            real :: sound_speed_avg, face_normal_speeds
            integer :: i , j, k
            
            do k = 1, kmx - 1
             do j = 1, jmx
              do i = 1, imx - 1
                sound_speed_avg = 0.5 * (sqrt(gm * y_pressure_left(i, j, k) / &
                                            y_density_left(i, j, k) ) + &
                                          sqrt(gm * y_pressure_right(i, j, k) / &
                                            y_density_right(i, j, k) ) )
                
                ! Compute '+' direction quantities
                face_normal_speeds = y_x_speed_left(i, j, k) * ynx(i, j, k) + &
                                     y_y_speed_left(i, j, k) * yny(i, j, k) + &
                                     y_z_speed_left(i, j, k) * ynz(i, j, k)
                M_perp_left = face_normal_speeds / sound_speed_avg
                alpha_plus = 0.5 * (1.0 + sign(1.0, M_perp_left))
                beta_left = -max(0, 1 - floor(abs(M_perp_left)))
                M_plus = 0.25 * ((1. + M_perp_left) ** 2.)
                D_plus = 0.25 * ((1. + M_perp_left) ** 2.) * (2. - M_perp_left)
                c_plus = (alpha_plus * (1.0 + beta_left) * M_perp_left) - &
                          beta_left * M_plus
                scrD_plus = (alpha_plus * (1. + beta_left)) - &
                        (beta_left * D_plus)

                ! First construct the F mass flux
                G_plus(1) = y_density_left(i, j, k) * sound_speed_avg * c_plus
                
                ! Is convective flux zero anywhere?
                if ((j .eq. 1) .or. (j .eq. jmx)) then
                    G_plus(1) = 0.
                end if
                ! Construct other fluxes in terms of the F mass flux
                G_plus(2) = (G_plus(1) * y_x_speed_left(i, j, k)) + &
                            (scrD_plus * y_pressure_left(i, j, k) * ynx(i, j, k))
                G_plus(3) = (G_plus(1) * y_y_speed_left(i, j, k)) + &
                            (scrD_plus * y_pressure_left(i, j, k) * yny(i, j, k))
                G_plus(4) = (G_plus(1) * y_z_speed_left(i, j, k)) + &
                            (scrD_plus * y_pressure_left(i, j, k) * ynz(i, j, k))
                G_plus(5) = G_plus(1) * &
                            ((0.5 * (y_x_speed_left(i, j, k) ** 2. + &
                                     y_y_speed_left(i, j, k) ** 2. + &
                                     y_z_speed_left(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * y_pressure_left(i, j, k) / &
                             y_density_left(i, j, k)))
         
                ! Multiply in the face areas
                G_plus(1) = G_plus(1) * yA(i, j, k)
                G_plus(2) = G_plus(2) * yA(i, j, k)
                G_plus(3) = G_plus(3) * yA(i, j, k)
                G_plus(4) = G_plus(4) * yA(i, j, k)
                G_plus(5) = G_plus(5) * yA(i, j, k)

                ! Compute '-' direction quantities

                face_normal_speeds = y_x_speed_right(i, j, k) * ynx(i, j, k) + &
                                     y_y_speed_right(i, j, k) * yny(i, j, k) + &
                                     y_z_speed_right(i, j, k) * ynz(i, j, k)
                M_perp_right = face_normal_speeds / sound_speed_avg
                alpha_minus = 0.5 * (1.0 - sign(1.0, M_perp_right))
                beta_right = -max(0, 1 - floor(abs(M_perp_right)))
                M_minus = - 0.25 * ((1. - M_perp_right) ** 2.)
                D_minus = 0.25 * ((1. - M_perp_right) ** 2.) * (2. + M_perp_right)
                c_minus = (alpha_minus * (1.0 + beta_right) * M_perp_right) - &
                          beta_right * M_minus
                scrD_minus = (alpha_minus * (1. + beta_right)) - &
                             (beta_right * D_minus)

                ! First construct the G mass flux
                G_minus(1) = y_density_right(i, j, k) * sound_speed_avg * c_minus
                
                ! Is convective flux zero anywhere?
                if ((j .eq. 1) .or. (j .eq. jmx)) then
                    G_minus(1) = 0.
                end if

                ! Construct other fluxes in terms of the G mass flux
                G_minus(2) = (G_minus(1) * y_x_speed_right(i, j, k)) + &
                             (scrD_minus * y_pressure_right(i, j, k) * ynx(i, j, k))
                G_minus(3) = (G_minus(1) * y_y_speed_right(i, j, k)) + &
                             (scrD_minus * y_pressure_right(i, j, k) * yny(i, j, k))
                G_minus(4) = (G_minus(1) * y_z_speed_right(i, j, k)) + &
                             (scrD_minus * y_pressure_right(i, j, k) * ynz(i, j, k))
                G_minus(5) = G_minus(1) * &
                             ((0.5 * (y_x_speed_right(i, j, k) ** 2. + &
                                      y_y_speed_right(i, j, k) ** 2. + &
                                      y_z_speed_right(i, j, k) ** 2.)) + &
                             ((gm / (gm - 1.)) * y_pressure_right(i, j, k) / &
                               y_density_right(i, j, k)))
         
                ! Multiply in the face areas
                G_minus(1) = G_minus(1) * yA(i, j, k)
                G_minus(2) = G_minus(2) * yA(i, j, k)
                G_minus(3) = G_minus(3) * yA(i, j, k)
                G_minus(4) = G_minus(4) * yA(i, j, k)
                G_minus(5) = G_minus(5) * yA(i, j, k)

                ! Get the total flux for a face
                G(i, j, k, :) = G_plus(:) + G_minus(:)
              end do
             end do
            end do 

          ! print *, 'Checking if G_conv flux is zero'
          ! print *, G(4, jmx, 2, 1)
          ! print *, G(20, 1, 1, 1)

        end subroutine compute_G_flux

        subroutine compute_H_flux()
            !-----------------------------------------------------------
            ! Compute the flux 'G' along the zeta faces
            !-----------------------------------------------------------
            
            implicit none

            real, dimension(1:n_var) :: H_plus, H_minus
            real :: M_perp_left, M_perp_right
            real :: alpha_plus, alpha_minus
            real :: beta_left, beta_right
            real :: M_plus, M_minus
            real :: D_plus, D_minus
            real :: c_plus, c_minus
            real :: scrD_plus, scrD_minus
            real :: sound_speed_avg, face_normal_speeds
            integer :: i , j, k
            
            do k = 1, kmx
             do j = 1, jmx - 1
              do i = 1, imx - 1
                sound_speed_avg = 0.5 * (sqrt(gm * z_pressure_left(i, j, k) / &
                                            z_density_left(i, j, k) ) + &
                                          sqrt(gm * z_pressure_right(i, j, k) / &
                                            z_density_right(i, j, k) ) )
                
                ! Compute '+' direction quantities
                face_normal_speeds = z_x_speed_left(i, j, k) * znx(i, j, k) + &
                                     z_y_speed_left(i, j, k) * zny(i, j, k) + &
                                     z_z_speed_left(i, j, k) * znz(i, j, k)
                M_perp_left = face_normal_speeds / sound_speed_avg
                alpha_plus = 0.5 * (1.0 + sign(1.0, M_perp_left))
                beta_left = -max(0, 1 - floor(abs(M_perp_left)))
                M_plus = 0.25 * ((1. + M_perp_left) ** 2.)
                D_plus = 0.25 * ((1. + M_perp_left) ** 2.) * (2. - M_perp_left)
                c_plus = (alpha_plus * (1.0 + beta_left) * M_perp_left) - &
                          beta_left * M_plus
                scrD_plus = (alpha_plus * (1. + beta_left)) - &
                        (beta_left * D_plus)

                ! First construct the H mass flux
                H_plus(1) = z_density_left(i, j, k) * sound_speed_avg * c_plus
                
                ! Is convective flux zero anywhere?
                if ((k .eq. 1) .or. (k .eq. kmx)) then
                    H_plus(1) = 0.
                end if

                ! Construct other fluxes in terms of the H mass flux
                H_plus(2) = (H_plus(1) * z_x_speed_left(i, j, k)) + &
                            (scrD_plus * z_pressure_left(i, j, k) * znx(i, j, k))
                H_plus(3) = (H_plus(1) * z_y_speed_left(i, j, k)) + &
                            (scrD_plus * z_pressure_left(i, j, k) * zny(i, j, k))
                H_plus(4) = (H_plus(1) * z_z_speed_left(i, j, k)) + &
                            (scrD_plus * z_pressure_left(i, j, k) * znz(i, j, k))
                H_plus(5) = H_plus(1) * &
                            ((0.5 * (z_x_speed_left(i, j, k) ** 2. + & 
                                     z_y_speed_left(i, j, k) ** 2. + &
                                     z_z_speed_left(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * z_pressure_left(i, j, k) / &
                             z_density_left(i, j, k)))
         
                ! Multiply in the face areas
                H_plus(1) = H_plus(1) * zA(i, j, k)
                H_plus(2) = H_plus(2) * zA(i, j, k)
                H_plus(3) = H_plus(3) * zA(i, j, k)
                H_plus(4) = H_plus(4) * zA(i, j, k)
                H_plus(5) = H_plus(5) * zA(i, j, k)

                ! Compute '-' direction quantities

                face_normal_speeds = z_x_speed_right(i, j, k) * znx(i, j, k) + &
                                     z_y_speed_right(i, j, k) * zny(i, j, k) + &
                                     z_z_speed_right(i, j, k) * znz(i, j, k)
                M_perp_right = face_normal_speeds / sound_speed_avg
                alpha_minus = 0.5 * (1.0 - sign(1.0, M_perp_right))
                beta_right = -max(0, 1 - floor(abs(M_perp_right)))
                M_minus = - 0.25 * ((1. - M_perp_right) ** 2.)
                D_minus = 0.25 * ((1. - M_perp_right) ** 2.) * (2. + M_perp_right)
                c_minus = (alpha_minus * (1.0 + beta_right) * M_perp_right) - &
                          beta_right * M_minus
                scrD_minus = (alpha_minus * (1. + beta_right)) - &
                             (beta_right * D_minus)

                ! First construct the F mass flux
                H_minus(1) = z_density_right(i, j, k) * sound_speed_avg * c_minus
                
                ! Is convective flux zero anywhere?
                if ((k .eq. 1) .or. (k .eq. kmx)) then
                    H_minus(1) = 0.
                end if

                ! Construct other fluxes in terms of the F mass flux
                H_minus(2) = (H_minus(1) * z_x_speed_right(i, j, k)) + &
                             (scrD_minus * z_pressure_right(i, j, k) * znx(i, j, k))
                H_minus(3) = (H_minus(1) * z_y_speed_right(i, j, k)) + &
                             (scrD_minus * z_pressure_right(i, j, k) * zny(i, j, k))
                H_minus(4) = (H_minus(1) * z_z_speed_right(i, j, k)) + &
                             (scrD_minus * z_pressure_right(i, j, k) * znz(i, j, k))
                H_minus(5) = H_minus(1) * &
                            ((0.5 * (z_x_speed_right(i, j, k) ** 2. + &
                                     z_y_speed_right(i, j, k) ** 2. + &
                                     z_z_speed_right(i, j, k) ** 2.)) + &
                            ((gm / (gm - 1.)) * z_pressure_right(i, j, k) / &
                             z_density_right(i, j, k)))
         
                ! Multiply in the face areas
                H_minus(1) = H_minus(1) * zA(i, j, k)
                H_minus(2) = H_minus(2) * zA(i, j, k)
                H_minus(3) = H_minus(3) * zA(i, j, k)
                H_minus(4) = H_minus(4) * zA(i, j, k)
                H_minus(5) = H_minus(5) * zA(i, j, k)

                ! Get the total flux for a face
                H(i, j, k, :) = H_plus(:) + H_minus(:)
              end do
             end do
            end do 
            
        end subroutine compute_H_flux

      ! subroutine compute_xi_face_quantities()
      !     !-----------------------------------------------------------
      !     ! Compute xi direction quantities used in F flux computation
      !     !
      !     ! x_c_plus and x_scrD_plus are used to calculate F_plus. 
      !     ! Similarly for F_minus. This subroutine computes these.
      !     !-----------------------------------------------------------

      !     implicit none
      !     real, dimension(imx, jmx-1) :: x_face_normal_speeds
      !     real, dimension(imx, jmx-1) :: M_plus, M_minus
      !     real, dimension(imx, jmx-1) :: D_plus, D_minus

      !     call dmsg(1, 'van_leer', 'compute_xi_face_quantities')

      !     x_sound_speed_avg = 0.5 * &
      !             (x_sound_speed_left() + x_sound_speed_right())

      !     ! Compute the '+' direction quantities
      !     x_face_normal_speeds = x_x_speed_left * xnx + x_y_speed_left * xny
      !     x_M_perp_left = x_face_normal_speeds / x_sound_speed_avg
      !     x_alpha_plus = 0.5 * (1.0 + sign(1.0, x_M_perp_left))
      !     x_beta_left = -max(0, 1 - floor(abs(x_M_perp_left)))
      !     M_plus = 0.25 * ((1. + x_M_perp_left) ** 2.)
      !     D_plus = 0.25 * ((1. + x_M_perp_left) ** 2.) * (2. - x_M_perp_left)
      !     x_c_plus = (x_alpha_plus * (1.0 + x_beta_left) * x_M_perp_left) - &
      !             x_beta_left * M_plus
      !     x_scrD_plus = (x_alpha_plus * (1. + x_beta_left)) - &
      !             (x_beta_left * D_plus)

      !     ! Compute the '-' direction quantities
      !     x_face_normal_speeds = x_x_speed_right * xnx + &
      !             x_y_speed_right * xny
      !     x_M_perp_right = x_face_normal_speeds / x_sound_speed_avg
      !     x_alpha_minus = 0.5 * (1.0 - sign(1.0, x_M_perp_right))
      !     x_beta_right = -max(0, 1 - floor(abs(x_M_perp_right)))
      !     M_minus = - 0.25 * ((1 - x_M_perp_right) ** 2.)
      !     D_minus = 0.25 * ((1 - x_M_perp_right) ** 2.) * &
      !             (2. + x_M_perp_right)
      !     x_c_minus = (x_alpha_minus * (1.0 + x_beta_right) * &
      !             x_M_perp_right) - (x_beta_right * M_minus)
      !     x_scrD_minus = (x_alpha_minus * (1.0 + x_beta_right)) - &
      !             (x_beta_right * D_minus)

      !     call dmsg(1, 'van_leer', 'compute_xi_face_quantities', 'Ended')

      ! end subroutine compute_xi_face_quantities

      ! function compute_F_plus() result(F_plus)

      !     implicit none
      !     real, dimension(imx, jmx-1, n_var) :: F_plus

      !     ! First construct the F mass flux
      !     F_plus(:, :, 1) = x_density_left * x_sound_speed_avg * x_c_plus
      !     ! Make the convective F mass flux zero at the right wall
      !     F_plus(imx, 32:63, 1) = 0
      !     
      !     ! Construct other fluxes in terms of the F mass flux
      !     F_plus(:, :, 2) = (F_plus(:, :, 1) * x_x_speed_left) + &
      !             (x_scrD_plus * x_pressure_left * xnx)
      !     F_plus(:, :, 3) = (F_plus(:, :, 1) * x_y_speed_left) + &
      !             (x_scrD_plus * x_pressure_left * xny)
      !     F_plus(:, :, 4) = F_plus(:, :, 1) * &
      !             ((0.5 * (x_x_speed_left ** 2. + x_y_speed_left ** 2.)) + &
      !             ((gm / (gm - 1.)) * x_pressure_left / x_density_left))

      !     ! Multiply in the face areas
      !     F_plus(:, :, 1) = F_plus(:, :, 1) * xA
      !     F_plus(:, :, 2) = F_plus(:, :, 2) * xA
      !     F_plus(:, :, 3) = F_plus(:, :, 3) * xA
      !     F_plus(:, :, 4) = F_plus(:, :, 4) * xA

      ! end function compute_F_plus

      ! function compute_F_minus() result(F_minus)
      !     
      !     implicit none
      !     real, dimension(imx, jmx-1, n_var) :: F_minus

      !     ! First construct the F mass flux
      !     F_minus(:, :, 1) = x_density_right * x_sound_speed_avg * x_c_minus
      !     ! Make the convective F mass flux zero at the right wall
      !     F_minus(imx, 32:63, 1) = 0
      !     
      !     ! Construct other fluxes in terms of the F mass flux
      !     F_minus(:, :, 2) = (F_minus(:, :, 1) * x_x_speed_right) + &
      !             (x_scrD_minus * x_pressure_right * xnx)
      !     F_minus(:, :, 3) = (F_minus(:, :, 1) * x_y_speed_right) + &
      !             (x_scrD_minus * x_pressure_right * xny)
      !     F_minus(:, :, 4) = F_minus(:, :, 1) * &
      !             ((0.5 * (x_x_speed_right ** 2. + x_y_speed_right ** 2.)) + &
      !             ((gm / (gm - 1.)) * x_pressure_right / x_density_right))

      !     ! Multiply in the face areas
      !     F_minus(:, :, 1) = F_minus(:, :, 1) * xA
      !     F_minus(:, :, 2) = F_minus(:, :, 2) * xA
      !     F_minus(:, :, 3) = F_minus(:, :, 3) * xA
      !     F_minus(:, :, 4) = F_minus(:, :, 4) * xA

      ! end function compute_F_minus

      ! subroutine compute_eta_face_quantities()
      !     !-----------------------------------------------------------
      !     ! Compute eta direction quantities used in G flux computation
      !     !
      !     ! y_c_plus and y_scrD_plus are used to calculate G_plus. 
      !     ! Similarly for G_minus. This subroutine computes these.
      !     !-----------------------------------------------------------

      !     implicit none
      !     real, dimension(imx-1, jmx) :: y_face_normal_speeds
      !     real, dimension(imx-1, jmx) :: M_plus, M_minus
      !     real, dimension(imx-1, jmx) :: D_plus, D_minus

      !     call dmsg(1, 'van_leer', 'compute_eta_face_quantities')

      !     y_sound_speed_avg = 0.5 * &
      !             (y_sound_speed_left() + y_sound_speed_right())

      !     ! Compute the '+' direction quantities
      !     y_face_normal_speeds = y_x_speed_left * ynx + y_y_speed_left * yny
      !     y_M_perp_left = y_face_normal_speeds / y_sound_speed_avg
      !     y_alpha_plus = 0.5 * (1.0 + sign(1.0, y_M_perp_left))
      !     y_beta_left = -max(0, 1 - floor(abs(y_M_perp_left)))
      !     M_plus = 0.25 * ((1. + y_M_perp_left) ** 2.)
      !     D_plus = 0.25 * ((1. + y_M_perp_left) ** 2.) * (2. - y_M_perp_left)
      !     y_c_plus = (y_alpha_plus * (1.0 + y_beta_left) * y_M_perp_left) - &
      !             y_beta_left * M_plus
      !     y_scrD_plus = (y_alpha_plus * (1. + y_beta_left)) - &
      !             (y_beta_left * D_plus)

      !     ! Compute the '-' direction quantities
      !     y_face_normal_speeds = y_x_speed_right * ynx + &
      !             y_y_speed_right(:, :) * yny
      !     y_M_perp_right = y_face_normal_speeds / y_sound_speed_avg
      !     y_alpha_minus = 0.5 * (1.0 - sign(1.0, y_M_perp_right))
      !     y_beta_right = -max(0, 1 - floor(abs(y_M_perp_right)))
      !     M_minus = - 0.25 * ((1 - y_M_perp_right) ** 2.)
      !     D_minus = 0.25 * ((1 - y_M_perp_right) ** 2.) * &
      !             (2. + y_M_perp_right)
      !     y_c_minus = (y_alpha_minus * (1.0 + y_beta_right) * &
      !             y_M_perp_right) - (y_beta_right * M_minus)
      !     y_scrD_minus = (y_alpha_minus * (1.0 + y_beta_right)) - &
      !             (y_beta_right * D_minus)
      !     call dmsg(1, 'van_leer', 'compute_eta_face_quantities', 'Ended')

      ! end subroutine compute_eta_face_quantities

      ! function compute_G_plus() result(G_plus)

      !     implicit none
      !     real, dimension(imx-1, jmx, 4) :: G_plus

      !     ! First construct the G mass flux
      !     G_plus(:, :, 1) = y_density_left * y_sound_speed_avg * y_c_plus
      !     ! Make the convective G mass flux zero at the wall
      !     G_plus(:, 1, 1) = 0
      !     G_plus(:, jmx, 1) = 0
      !     ! Construct other fluxes in terms of the G mass flux
      !     G_plus(:, :, 2) = (G_plus(:, :, 1) * y_x_speed_left) + &
      !             (y_scrD_plus * y_pressure_left * ynx)
      !     G_plus(:, :, 3) = (G_plus(:, :, 1) * y_y_speed_left) + &
      !             (y_scrD_plus * y_pressure_left * yny)
      !     G_plus(:, :, 4) = G_plus(:, :, 1) * &
      !             ((0.5 * (y_x_speed_left ** 2. + y_y_speed_left ** 2.)) + &
      !             ((gm / (gm - 1.)) * y_pressure_left / y_density_left))

      !     ! Multiply in the face areas
      !     G_plus(:, :, 1) = G_plus(:, :, 1) * yA
      !     G_plus(:, :, 2) = G_plus(:, :, 2) * yA
      !     G_plus(:, :, 3) = G_plus(:, :, 3) * yA
      !     G_plus(:, :, 4) = G_plus(:, :, 4) * yA

      ! end function compute_G_plus

      ! function compute_G_minus() result(G_minus)
      !     
      !     implicit none
      !     real, dimension(imx-1, jmx, 4) :: G_minus

      !     ! First construct the G mass flux
      !     G_minus(:, :, 1) = y_density_right * y_sound_speed_avg * y_c_minus
      !     ! Make the convective G mass flux zero at the wall
      !     G_minus(:, 1, 1) = 0
      !     G_minus(:, jmx, 1) = 0
      !     ! Construct other fluxes in terms of the G mass flux
      !     G_minus(:, :, 2) = (G_minus(:, :, 1) * y_x_speed_right) + &
      !             (y_scrD_minus * y_pressure_right * ynx)
      !     G_minus(:, :, 3) = (G_minus(:, :, 1) * y_y_speed_right) + &
      !             (y_scrD_minus * y_pressure_right * yny)
      !     G_minus(:, :, 4) = G_minus(:, :, 1) * &
      !             ((0.5 * (y_x_speed_right ** 2. + &
      !                 y_y_speed_right ** 2.)) + &
      !             ((gm / (gm - 1.)) * y_pressure_right / y_density_right))

      !     ! Multiply in the face areas
      !     G_minus(:, :, 1) = G_minus(:, :, 1) * yA
      !     G_minus(:, :, 2) = G_minus(:, :, 2) * yA
      !     G_minus(:, :, 3) = G_minus(:, :, 3) * yA
      !     G_minus(:, :, 4) = G_minus(:, :, 4) * yA

      ! end function compute_G_minus

      ! function compute_residue() result(residue)
      !     !-----------------------------------------------------------
      !     ! Compute the residue using the Van-Leer scheme
      !     !-----------------------------------------------------------

      !     implicit none
      !     real, dimension(1:imx-1, 1:jmx-1, 1:kmx-1, 1:n_var) :: residue
      !     integer :: i, j, k, l

      !     call dmsg(1, 'van_leer', 'compute_residue')

      !     call compute_F_flux()
      !     if (any(isnan(F))) then
      !         call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: F flux Nan detected')
      !         stop
      !     end if    
      !     
      !     call compute_G_flux()
      !     if (any(isnan(G))) then 
      !         call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: G flux Nan detected')
      !         stop
      !     end if    
      !     
      !     call compute_H_flux()
      !     if (any(isnan(H))) then
      !         call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: H flux Nan detected')
      !         stop
      !     end if
      !     
      !     do l = 1, n_var
      !      do k = 1, kmx - 1
      !       do j = 1, jmx - 1
      !        do i = 1, imx - 1
      !        residue(i, j, k, l) = F(i+1, j, k, l) - F(i, j, k, l) &
      !                            + G(i, j+1, k, l) - G(i, j, k, l) &
      !                            + H(i, j, k+1, l) - H(i, j, k, l)
      !        end do
      !       end do
      !      end do
      !     end do
      ! 
      ! end function compute_residue

      ! function get_residue() result(residue)
      !     !-----------------------------------------------------------
      !     ! Return the VL residue
      !     !-----------------------------------------------------------
      !     
      !     implicit none
      !     real, dimension(imx-1, jmx-1, kmx-1, n_var) :: residue

      !     residue = compute_residue()

      ! end function get_residue

        subroutine compute_fluxes()
            
            implicit none
            
            call dmsg(1, 'van_leer', 'compute_fluxes')

            call compute_F_flux()
            if (any(isnan(F))) then
                call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: F flux Nan detected')
                stop
            end if    

            call compute_G_flux()
            if (any(isnan(G))) then 
                call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: G flux Nan detected')
                stop
            end if    
            
            call compute_H_flux()
            if (any(isnan(H))) then
                call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: H flux Nan detected')
                stop
            end if

        end subroutine compute_fluxes

        function get_residue() result(residue)
            !-----------------------------------------------------------
            ! Compute the residue using the Van-Leer scheme
            !-----------------------------------------------------------

            implicit none
            
            integer :: i, j, k, l
            real, dimension(1:imx-1, 1:jmx-1, 1:kmx-1, 1:n_var) :: residue

            call dmsg(1, 'van_leer', 'compute_residue')

 !          call compute_F_flux()
 !          if (any(isnan(F))) then
 !              call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: F flux Nan detected')
 !              stop
 !          end if    
 !          
 !          call compute_G_flux()
 !          if (any(isnan(G))) then 
 !              call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: G flux Nan detected')
 !              stop
 !          end if    
 !          
 !          call compute_H_flux()
 !          if (any(isnan(H))) then
 !              call dmsg(5, 'van_leer', 'compute_residue', 'ERROR: H flux Nan detected')
 !              stop
 !          end if
            
            do l = 1, n_var
             do k = 1, kmx - 1
              do j = 1, jmx - 1
               do i = 1, imx - 1
               residue(i, j, k, l) = F(i+1, j, k, l) - F(i, j, k, l) &
                                   + G(i, j+1, k, l) - G(i, j, k, l) &
                                   + H(i, j, k+1, l) - H(i, j, k, l)
               end do
              end do
             end do
            end do
        
        end function get_residue

end module van_leer
