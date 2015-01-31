module physics
    ! always keep module named 'physics'
    !   even if the file name changes.
    !   - this lets you create different physics for different runs,
    !       and have things like simply by compiling the right physics file
    !       !!! NOTE !!! 
    !           - this means you MUST make sure you have a recent,
    !             compiled version of the right physics.
    
    implicit none
    private calc_energy, calc_angular_momentum_scalar
    public calc_force, accuracy_check

    type, public :: accuracy_state
            real(kind=16) :: energy, angular_momentum
        contains
            procedure, pass :: initialize => accuracy_state_initialize
    end type accuracy_state

    real(kind=16), parameter     :: energy_tol           = 1d-11
    real(kind=16), parameter     :: angular_momentum_tol = 1d-6

    contains

subroutine calc_force(n_bodies, x, mass, force)
    ! Calculate gravitational forces given current positions

    ! inputs:
    !   n_bodies    - integer - number of interacting bodies
    !   x           - (n_bodies x 3) float array - positions 

    ! outputs:
    !   f           - (n_bodies x 3) float array - forces 

    ! side effects:
    !   none (this is enforced)

    ! assumes:

    ! features for later:

    implicit none
    integer,                               intent(in)  :: n_bodies
    real(kind=16), dimension(n_bodies, 3), intent(in)  :: x
    real(kind=16), dimension(n_bodies   ), intent(in)  :: mass
    real(kind=16), dimension(n_bodies, 3), intent(out) :: force
    real(kind=16), dimension(3)                        :: df
    integer                                            :: i,j ! loop variables

    force = 0.

    do concurrent (i=1:n_bodies, j=1:n_bodies, i.gt.j )
        df = mass(i)*mass(j) &
            * ( x(j,:) - x(i,:) ) &   ! displacement vector from j to i
            / sqrt(sum( (x(j,:) - x(i,:))**2 ))**3 ! distance^3
        force(i,:) = force(i,:) + df
        force(j,:) = force(j,:) - df ! "-" sign for flip of displacement
    end do

end subroutine calc_force
    pure function calc_energy(n_bodies, x, v, mass) result(energy)
        implicit none
        integer,                               intent(in)  :: n_bodies
        real(kind=16), dimension(n_bodies, 3), intent(in)  :: x, v
        real(kind=16), dimension(n_bodies   ), intent(in)  :: mass
        real(kind=16)                                      :: energy
        integer                                            :: i,j ! loop variables

        energy = 0.

        do concurrent (i=1:n_bodies, j=1:3)
            energy = energy + 1./2 * mass(i) * v(i,j)**2
        end do

        do concurrent (i=1:n_bodies, j=1:n_bodies, i.gt.j )
            energy  = energy -  mass(i)*mass(j) &
                / sqrt(sum((/ x(j,:) - x(i,:) /)**2))   ! distance
        end do
        
    end function calc_energy

    pure function calc_angular_momentum_scalar(n_bodies, x, v, mass) result(angular_momentum)
        ! only calculates angular momentum in z-hat direction
        implicit none
        integer,                               intent(in)  :: n_bodies
        real(kind=16), dimension(n_bodies, 3), intent(in)  :: x, v
        real(kind=16), dimension(n_bodies   ), intent(in)  :: mass
        real(kind=16)                                      :: angular_momentum
        integer                                            :: i ! loop variables

        angular_momentum = 0.

        do concurrent (i=1:n_bodies)
            ! ! ! L = m * (r x v) = m * (x*v_y - y*v_x)
            angular_momentum = angular_momentum &
                + mass(i) * ( (x(i,1) * v(i,2)) - (x(i,2) * v(i,1)) )
        end do

    end function calc_angular_momentum_scalar

    function accuracy_check(state_initial, state_current) result (check)
        implicit none
        class(accuracy_state), intent(in) :: state_initial, state_current
        logical :: check
        real(kind=16) :: energy_err, angular_momentum_err

        energy_err = abs( (state_current%energy - state_initial%energy) &
            / (state_initial%energy) )

        angular_momentum_err = abs( state_current%angular_momentum &
            - state_initial%angular_momentum )

        check = .False.
        if (energy_err.lt.energy_tol) then
            if (angular_momentum_err.lt.angular_momentum_tol) then
                check = .True.
            end if
        end if

    end function

    subroutine accuracy_state_initialize(state, n_bodies, x, v, mass)
        implicit none
        class(accuracy_state),                 intent(out) :: state
        integer,                               intent(in)  :: n_bodies
        real(kind=16), dimension(n_bodies, 3), intent(in)  :: x, v
        real(kind=16), dimension(n_bodies   ), intent(in)  :: mass
        real(kind=16) :: energy, angular_momentum

        energy = calc_energy(n_bodies, x, v, mass)
        angular_momentum = calc_angular_momentum_scalar(n_bodies, x, v, mass)

        state%energy = energy
        state%angular_momentum = angular_momentum

    end subroutine


    end module