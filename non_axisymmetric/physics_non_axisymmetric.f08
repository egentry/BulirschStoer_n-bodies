module physics
    ! always keep module named 'physics'
    !   even if the file name changes.
    !   - this lets you create different physics for different runs,
    !       and have things like simply by compiling the right physics file
    !       !!! NOTE !!! 
    !           - this means you MUST make sure you have a recent,
    !             compiled version of the right physics.

    implicit none
    private calc_energy
    public calc_force, accuracy_check
    integer, parameter           :: precision=8 ! double precision

    type, public :: accuracy_state
            real(kind=precision) :: energy
        contains
            procedure, pass :: initialize => accuracy_state_initialize
    end type accuracy_state


    real(kind=precision), parameter     :: energy_tol           = 1d-4

    real(kind=precision), parameter     :: eccentricity    = .3
    real(kind=precision), parameter     :: energy_initial  = -1
    real(kind=precision), parameter     :: x_max           = exp(energy_initial) &
                                                        / (1 - eccentricity)



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
        integer,                               intent(in)         :: n_bodies
        real(kind=precision), dimension(n_bodies, 3), intent(in)  :: x
        real(kind=precision), dimension(n_bodies   ), intent(in)  :: mass
        real(kind=precision), dimension(n_bodies, 3), intent(out) :: force
        real(kind=precision)                                      :: denominator
        real(kind=precision)                                      :: x1, x2 ! x,y coords
        integer                                                   :: i,j ! loop variables

        do concurrent (i=1:n_bodies)
            x1 = x(i,1)
            x2 = x(i,2)

            denominator = ( (x1**2 + x2**2) &
                - eccentricity * x1 * sqrt(x1**2 + x2**2) )

            force(i,1) = -1 * mass(i) * ( (x1 - eccentricity * sqrt(x1**2 + x2**2)) &
                / denominator )
            force(i,2) = -1 * mass(i) * ( (x2) &
                / denominator )
            force(i,3) = 0.
        end do

    end subroutine calc_force


    pure function calc_energy(n_bodies, x, v, mass) result(energy)
        implicit none
        integer,                               intent(in)         :: n_bodies
        real(kind=precision), dimension(n_bodies, 3), intent(in)  :: x, v
        real(kind=precision), dimension(n_bodies   ), intent(in)  :: mass
        real(kind=precision)                                      :: energy
        integer                                                   :: i,j ! loop variables

        energy = 0.

        do concurrent (i=1:n_bodies, j=1:3)
            energy = energy + 1./2 * mass(i) * v(i,j)**2
        end do

        do concurrent (i=1:n_bodies)
            energy  = energy + log( sqrt(x(i,1)**2 + x(i,2)**2) &
                - eccentricity * x(i,1) )
        end do
        
    end function calc_energy

    subroutine calc_initial_velocity(n_bodies, x, v)
        ! Sets x velocity = 0, then:
        !
        ! Calculates the initial y velocities,
        !   assuming: 
        !       - 0 < x < x_max
        !           - where x_max is maximum allowed at energy_initial=-1
        !       - x velocity = 0
        !       - y          = 0
        !       - z          = 0 (always)
        !       - z velocity = 0

        ! inputs:
        !   n_bodies    - integer - number of interacting bodies
        !   x           - (n_bodies x 3) float array - positions 

        ! outputs:
        !              - (n_bodies x 3) float array - velocities 

        ! side effects:
        !   overwrites any value that existed in velocities

        ! assumes:
        !   - see above

        ! features for later:

        implicit none
        integer,                               intent(in)  :: n_bodies
        real(kind=precision), dimension(n_bodies, 3), intent(in)  :: x
        real(kind=precision), dimension(n_bodies, 3), intent(out) :: v
        real(kind=precision)                                      :: x0
        integer                                            :: i,j ! loop variables

        do concurrent (i=1:n_bodies)
            v(i,1) = 0.
            v(i,2) = sqrt( 2 * (energy_initial  &
                - log( x(i,1) * (1 - eccentricity) )) )
            v(i,3) = 0.
        end do

    end subroutine calc_initial_velocity

    function accuracy_check(state_initial, state_current) result (check)
        implicit none
        class(accuracy_state), intent(in) :: state_initial, state_current
        logical :: check
        real(kind=precision) :: energy_err

        energy_err = abs( (state_current%energy - state_initial%energy) &
            / (state_initial%energy) )

        check = .False.
        if (energy_err.lt.energy_tol) then
            check = .True.
        end if

    end function

    subroutine accuracy_state_initialize(state, n_bodies, x, v, mass)
        implicit none
        class(accuracy_state),                 intent(out) :: state
        integer,                               intent(in)  :: n_bodies
        real(kind=precision), dimension(n_bodies, 3), intent(in)  :: x, v
        real(kind=precision), dimension(n_bodies   ), intent(in)  :: mass
        real(kind=precision) :: energy

        energy = calc_energy(n_bodies, x, v, mass)

        state%energy = energy

    end subroutine


end module