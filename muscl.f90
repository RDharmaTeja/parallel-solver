module muscl
    !-----------------------------------------------------------------
    ! MUSCL (Monotone Upwing Schemes for Scalar Conservation Laws is
    ! a scheme which replaces the piecewise constant approximation by
    ! reconstructing the states at the left and right side of each face.
    ! This is a one parameter upwind scheme which results in at most 3rd
    ! order accuracy.
    !
    ! The MUSCL scheme alone creates non-physical oscillations near 
    ! discontinuities like shocks. Hence, MUSCL is combined with
    ! some TVD (Total Variation Diminishing) to reduce such oscillations.
    ! TVD schemes also ensure that no new extrema of the state variables
    ! is created at the faces.
    !-----------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx, kmx
    use state, only: qp, n_var

    implicit none
    private

    ! Private variables
    real, dimension(:, :, :, :), allocatable, target :: x_qp_left, &
        x_qp_right, y_qp_left, y_qp_right, z_qp_left, z_qp_right
    real :: phi, kappa
!   character(len=30) :: TVD_scheme

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_muscl_states
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right
    public :: z_qp_left, z_qp_right

 !  TVD_scheme = trim('koren')

    contains
        
        subroutine setup_scheme()

        implicit none

        call dmsg(1, 'muscl', 'setup_muscl')

        phi = 1.0

        kappa = -1.0

        call alloc(x_qp_left, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_left.')
        call alloc(x_qp_right, 1, imx, 1, jmx-1, 1, kmx-1, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_right.')

        call alloc(y_qp_left, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_left.')
        call alloc(y_qp_right, 1, imx-1, 1, jmx, 1, kmx-1, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_right.')

        call alloc(z_qp_left, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_left.')
        call alloc(z_qp_right, 1, imx-1, 1, jmx-1, 1, kmx, 1, n_var, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_right.')

        end subroutine setup_scheme


        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'muscl', 'destroy_muscl')

            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)
            call dealloc(z_qp_left)
            call dealloc(z_qp_right)

        end subroutine destroy_scheme


!       function min_mod(r)
!           
!           implicit none

!           real, intent(in) :: r
!           real :: min_mod
!           min_mod = min(1., (3 - kappa) * r / (1 - kappa))
!       
!       end function min_mod

!       function koren(r)

!           implicit none
!           
!           real, intent(in) :: r
!           real :: koren
!           koren = max(0., min(2*r, (2 + r)/3., 2.))

!       end function koren(r)


        subroutine compute_xi_face_states()

            implicit none

            integer :: i, j, k, l
            real :: psi1, psi2, fd, bd, r

            phi = 1.0
            kappa = -1.0

            do l = 1, n_var
             do k = 1, kmx - 1
              do j = 1, jmx - 1
               do i = 1, imx - 1
                ! All faces interior only (even at boundaries)
                ! Hence (i=1, left) and (i=imx, right) will be dealt separately
                ! Koren limiter for now
                ! From paper: delta: forward difference 'fd'
                !             nabla: backward difference 'bd'
                !TODO: Generalise TVD scheme functions
                fd = qp(i+1, j, k, l) - qp(i, j, k, l)
                bd = qp(i, j, k, l) - qp(i-1, j, k, l)
                r = fd / bd
             !  psi1 = min(1., (3 - kappa) * r / (1 - kappa))
                psi1 = max(0., min(2*r, (2 + r)/3., 2.))
                r = bd / fd
             !  psi2 = min(1., (3 - kappa) * r / (1 - kappa))
                psi2 = max(0., min(2*r, (2 + r)/3., 2.))

                x_qp_left(i+1, j, k, l) = qp(i, j, k, l) + 0.25*phi* &
                    (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
                x_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                    (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
               end do
              end do
             end do
            end do

            do k = 1, kmx - 1
             do j = 1, jmx - 1
                ! Exterior boundaries
                x_qp_left(1, j, k, :) = 0.5 * (qp(0, j, k, :) + qp(1, j, k, :))
                x_qp_right(imx, j, k, :) = 0.5 * (qp(imx-1, j, k, :) + qp(imx, j, k, :))
             end do
            end do

        end subroutine compute_xi_face_states


        subroutine compute_eta_face_states()

            implicit none

            integer :: i, j, k, l
            real :: psi1, psi2, fd, bd, r

            do l = 1, n_var
             do k = 1, kmx - 1
              do j = 1, jmx - 1
               do i = 1, imx - 1
                ! All faces interior only (even at boundaries)
                ! Hence (j=1, left) and (j=jmx, right) will be dealt separately
                ! Koren limiter for now
                ! From paper: delta: forward difference 'fd'
                !             nabla: backward difference 'bd'
                !TODO: Generalise TVD scheme functions
                fd = qp(i, j+1, k, l) - qp(i, j, k, l)
                bd = qp(i, j, k, l) - qp(i, j-1, k, l)
                r = fd / bd
                psi1 = max(0., min(2*r, (2 + r)/3., 2.))
             !  psi1 = min(1., (3 - kappa) * r / (1 - kappa))
                r = bd / fd
                psi2 = max(0., min(2*r, (2 + r)/3., 2.))
             !  psi2 = min(1., (3 - kappa) * r / (1 - kappa))

                y_qp_left(i, j+1, k, l) = qp(i, j, k, l) + 0.25*phi* &
                    (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
                y_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                    (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
               end do
              end do
             end do
            end do

            do k = 1, kmx - 1
             do i = 1, imx - 1
                ! Exterior boundaries
                y_qp_left(i, 1, k, :) = 0.5 * (qp(i, 0, k, :) + qp(i, 1, k, :))
                y_qp_right(i, jmx, k, :) = 0.5 * (qp(i, jmx-1, k, :) + qp(i, jmx, k, :))
             end do
            end do

        end subroutine compute_eta_face_states


        subroutine compute_zeta_face_states()

            implicit none

            real :: psi1, psi2, fd, bd, r
            integer :: i, j, k, l

            do k = 1, kmx - 1
             do l = 1, n_var
              do j = 1, jmx - 1
               do i = 1, imx - 1
                ! All faces interior only (even at boundaries)
                ! Hence (k=1, left) and (k=kmx, right) will be dealt separately
                ! Koren limiter for now
                ! From paper: delta: forward difference 'fd'
                !             nabla: backward difference 'bd'
                !TODO: Generalise TVD scheme functions
                fd = qp(i, j, k+1, l) - qp(i, j, k, l)
                bd = qp(i, j, k, l) - qp(i, j, k-1, l)
                r = fd / bd
                psi1 = max(0., min(2*r, (2 + r)/3., 2.))
             !  psi1 = min(1., (3 - kappa) * r / (1 - kappa))
                r = bd / fd
                psi2 = max(0., min(2*r, (2 + r)/3., 2.))
             !  psi2 = min(1., (3 - kappa) * r / (1 - kappa))

                z_qp_left(i, j, k+1, l) = qp(i, j, k, l) + 0.25*phi* &
                    (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
                z_qp_right(i, j, k, l) = qp(i, j, k, l) - 0.25*phi* &
                    (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
               end do
              end do
             end do
            end do

            do j = 1, jmx - 1
             do i = 1, imx - 1
                ! Exterior boundaries
                z_qp_left(i, j, 1, :) = 0.5 * (qp(i, k, 0, :) + qp(i, j, 1, :))
                z_qp_right(i, j, kmx, :) = 0.5 * (qp(i, j, kmx-1, :) + qp(i, j, kmx, :))
             end do
            end do

        end subroutine compute_zeta_face_states


        subroutine compute_muscl_states()
            !---------------------------------------------------------
            ! Implement MUSCL scheme to get left and right states at
            ! each face
            !---------------------------------------------------------
            
            call compute_xi_face_states()
            call compute_eta_face_states()
            call compute_zeta_face_states()

        end subroutine compute_muscl_states

end module muscl
