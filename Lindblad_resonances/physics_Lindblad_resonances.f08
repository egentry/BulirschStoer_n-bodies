module physics
    ! always keep module named 'physics'
    !   even if the file name changes.
    !   - this lets you create different physics for different runs,
    !       and have things like simply by compiling the right physics file
    !       !!! NOTE !!! 
    !           - this means you MUST make sure you have a recent,
    !             compiled version of the right physics.

    implicit none
    public calc_energy
    public calc_force, accuracy_check
    integer, parameter           :: precision=8 ! double precision

    type, public :: accuracy_state
            real(kind=precision) :: energy
        contains
            procedure, pass :: initialize => accuracy_state_initialize
    end type accuracy_state

    logical, parameter                  :: allow_z_motion = .False. ! for bulirsch_stoer


    real(kind=precision), parameter     :: energy_tol      = 1d-4

    real(kind=precision), parameter     :: epsilon         = .1
    real(kind=precision), parameter     :: R_c             = 0.5
    real(kind=precision), parameter     :: V_0             = 1
    real(kind=precision), parameter     :: Omega_b         = 1.25
    integer,              parameter     :: m               = 2


    ! ! ! DERIVED CONSTANTS
    ! location of resonances (inner Lindblad does not exist)
    real(kind=precision), parameter     :: R_corot         = 0.6245
    real(kind=precision), parameter     :: R_outer         = 1.308363879053123
    real(kind=precision), parameter     :: R_general       = 1

    ! unperturbed tangential velocities at resonances, in rotating frame
    real(kind=precision), parameter     :: v_corot_tang    =  0
    real(kind=precision), parameter     :: v_outer_tang    = -0.70134187239444
    real(kind=precision), parameter     :: v_general_tang  = -0.35557280900008




    contains

    subroutine calc_force(n_bodies, x, v, mass, force)
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
        integer,                                      intent(in)  :: n_bodies
        real(kind=precision), dimension(n_bodies, 3), intent(in)  :: x
        real(kind=precision), dimension(n_bodies, 3), intent(in)  :: v
        real(kind=precision), dimension(n_bodies   ), intent(in)  :: mass
        real(kind=precision), dimension(n_bodies, 3), intent(out) :: force
        real(kind=precision)                                      :: dLog_term
        real(kind=precision)                                      :: dcos_term
        real(kind=precision)                                      :: x1, x2 ! x,y coords
        real(kind=precision)                                      :: R_squared
        integer                                                   :: i ! loop variables

        do concurrent (i=1:n_bodies)
            x1 = x(i,1)  ! x coord
            x2 = x(i,2)  ! y coord
            R_squared = x1**2 + x2**2   ! radial coordinate, SQUARED

            dLog_term = V_0**2 * (1 + (epsilon**2 * ((x1**2 - x2**2) / R_squared)) ) &
                / (R_c**2 + R_squared)

            dcos_term = 2 * V_0**2 * log(R_c**2 + R_squared) &
                * (epsilon**2 * x1 * x2) / (R_squared**2)


            force(i,1) = mass(i) * ( -x1*dLog_term - x2*dcos_term )
            force(i,2) = mass(i) * ( -x2*dLog_term + x1*dcos_term )
            force(i,3) = 0.

            ! Add centrifugal:
            force(i,1) = force(i,1) + mass(i) * Omega_b**2 * x(i,1)
            force(i,2) = force(i,2) + mass(i) * Omega_b**2 * x(i,2)
            force(i,3) = force(i,3) + 0.          

!             ! Add coriolis:
            force(i,1) = force(i,1) + 2 * mass(i) * Omega_b * v(i,2)
            force(i,2) = force(i,2) - 2 * mass(i) * Omega_b * v(i,1)
            force(i,3) = force(i,3) + 0.                 
        end do

    end subroutine calc_force


    pure function calc_energy(n_bodies, x, v, mass) result(energy)
        implicit none
        integer,                                      intent(in)  :: n_bodies
        real(kind=precision), dimension(n_bodies, 3), intent(in)  :: x, v
        real(kind=precision), dimension(n_bodies   ), intent(in)  :: mass
        real(kind=precision)                                      :: energy
        integer                                                   :: i,j ! loop variables

        energy = 0.

        do concurrent (i=1:n_bodies, j=1:3)
            energy = energy + 1./2 * mass(i) * v(i,j)**2
        end do

        do concurrent (i=1:n_bodies)
            energy  = energy + mass(i) * (1./2) * V_0**2 &
                * log(R_c**2 + x(i,1)**2 + x(i,2)**2) &
                * (1 + epsilon**2  &
                    * ((x(i,1)**2 - x(i,2)**2) / (x(i,1)**2 + x(i,2)**2)) )
            energy = energy - ( mass(i)*Omega_b**2 * (x(i,1)**2 + x(i,2)**2)/2 )  ! centrifugal barrier
        end do
        
    end function calc_energy


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
        class(accuracy_state),                        intent(out) :: state
        integer,                                      intent(in)  :: n_bodies
        real(kind=precision), dimension(n_bodies, 3), intent(in)  :: x, v
        real(kind=precision), dimension(n_bodies   ), intent(in)  :: mass
        real(kind=precision) :: energy

        energy = calc_energy(n_bodies, x, v, mass)

        state%energy = energy

    end subroutine


end module