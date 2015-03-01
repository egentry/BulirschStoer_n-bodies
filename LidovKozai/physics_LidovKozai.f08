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

    logical, parameter              :: allow_z_motion = .True. ! for bulirsch_stoer
    real(kind=precision), parameter :: pi = 4 * atan(1.)

    
    ! ! !  TUNABLE QUANTITIES

    real(kind=precision), parameter :: energy_tol  = 1d-6

    real(kind=precision), parameter :: M_star      = 1
    real(kind=precision), parameter :: M_planet    = 0.01
    real(kind=precision), parameter :: M_outer     = 1

    ! initial semi-major axes from central star
    real(kind=precision), parameter :: a_planet_0  = 1
    real(kind=precision), parameter :: a_outer_0   = 100


    real(kind=precision), parameter :: Theta       = .25 ! conserved integral of motion
    real(kind=precision), parameter :: x_0         = .4


    ! ! ! DERIVED QUANTITIES
    real(kind=precision), parameter :: e_0         = sqrt(1 - x_0)
    real(kind=precision), parameter :: cosi_0      = sqrt(Theta / x_0) !cosine of initial inclination
    real(kind=precision), parameter :: i_0         = acos(cosi_0)
    
    ! initial distances of bodies, from central star
    real(kind=precision), parameter :: r_planet_0 = a_planet_0 * (1 - e_0)
    real(kind=precision), parameter :: r_outer_0  = a_outer_0

    ! periods of unperturbed orbits
    real(kind=precision), parameter :: P_planet   = 2*pi*sqrt(a_planet_0**3 / (M_star + M_planet))
    real(kind=precision), parameter :: P_outer    = 2*pi*sqrt(a_outer_0**3  / (M_star + M_outer ))

    !angular momentum magnitudes of unperturbed orbits
    real(kind=precision), parameter :: l_planet   = 2*pi*a_planet_0**2 * sqrt(1 - e_0**2) / P_planet
    real(kind=precision), parameter :: l_outer    = 2*pi*a_outer_0**2                     / P_outer
 
    ! initial velocities of orbits
    real(kind=precision), parameter :: v_planet_0 = l_planet / r_planet_0
    real(kind=precision), parameter :: v_outer_0  = l_outer  / r_outer_0


    contains

    subroutine calc_force(n_bodies, x, v, mass, force)
        ! Calculate gravitational forces given current positions

        ! inputs:
        !   n_bodies    - integer - number of interacting bodies
        !   x           - (n_bodies x 3) float array - positions 
        !   v           - (n_bodies x 3) float array - velocities 
        !               - not used

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
        real(kind=precision), dimension(3)                        :: df
        integer                                                   :: i,j ! loop variables

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
        integer,                                      intent(in)  :: n_bodies
        real(kind=precision), dimension(n_bodies, 3), intent(in)  :: x, v
        real(kind=precision), dimension(n_bodies   ), intent(in)  :: mass
        real(kind=precision)                                      :: energy
        integer                                                   :: i,j ! loop variables

        energy = 0.

        do concurrent (i=1:n_bodies, j=1:3)
            energy = energy + 1./2 * mass(i) * v(i,j)**2
        end do

        do concurrent (i=1:n_bodies, j=1:n_bodies, i.gt.j )
            energy  = energy -  mass(i)*mass(j) &
                / sqrt(sum((/ x(j,:) - x(i,:) /)**2))   ! distance
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