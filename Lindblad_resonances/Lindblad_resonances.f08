program Lindblad_resonances
    use physics
    use bulirsch_stoer
    implicit none

    ! the main things you need to set are:
    !    In physics_Lindblad.f08:
    !       - precision - sets float precision (e.g. 8 = double)
    !       - energy_tol - fractional tolerance for energy conservation *each step*

    !       - eccentricity - eccentricity of the potential well
    !       - R_c          - scale radius
    !       - V_0          - scale velocity
    !       - Omega_b      - rotational frequency of perturbation pattern
    !       - m            - number of azimuthal modes in perturbation
    !                      - m=2 for bar
    !       (R_corot, R_outer, v_corot_tang, v_outer_tang must also be changed
    !           to fit the given values in R_c, V_0, etc.
    !           See Binney & Tremaine Sec 3.3.3 to derive these values)

    !   In this file:
    !       - n_bodies
    !           - only two resonances exist for current potential
    !           - for more bodies, you'll need to hard-code initial conditions
    !       - t_max
    !           - how long should the integration continue?
    !       - t_print
    !           - how often should you print?


    integer, parameter                              :: n_bodies = 2
    ! 3-position, 3-velocity
    real(kind=precision), dimension(n_bodies, 3)    :: x, v
    real(kind=precision), dimension(n_bodies)       :: mass

    real(kind=precision), parameter     :: H_big = 1./real(10**(4), kind=precision) ! large step size
    real(kind=precision)                :: t
    real(kind=precision), parameter     :: t_max   = 5000.
    real(kind=precision), parameter     :: t_print = .001  ! print every t_print,

    ! defined in the physics module -- changes for each potential
    type(accuracy_state)                :: state_initial, state_final

    integer                             :: i,j,k        ! loop variables

    integer                             :: save_file_unit
    character(len=6)                    :: ith_body_str ! tmp variable used in filenames
    character(len=15)                   :: filename     ! tmp variable used in filenames
    character(len=15)                   :: data_format  ! for saving

    ! INITIALIZE

    ! initialize x_0 to a range of 0 < x_0 < x_max
    x(1,1) = R_corot
    x(1,2) = 0
    x(1,3) = 0

    x(2,1) = R_outer
    x(2,2) = 0
    x(2,3) = 0

    v(1,1) = v_corot_tang * epsilon
    v(1,2) = v_corot_tang
    v(1,3) = 0

    v(2,1) = v_outer_tang * epsilon
    v(2,2) = v_outer_tang
    v(2,3) = 0

    mass = 1.

    call state_initial%initialize(n_bodies, x, v, mass)

    print *, 'energy (initial): '
    print *,  state_initial%energy


    ! ! ! BEGIN INTEGRATIONS
    data_format = "(7e13.5)"
    do i = 1,n_bodies
        print *, 'starting body: '
        print *, i
        write (ith_body_str, "(I4.4)") i
        filename = 'body_'//trim(ith_body_str)//'.dat'

        open(newunit=save_file_unit, file=filename, status='unknown')
        write(save_file_unit, *) &
            "t",achar(9),achar(9),achar(9),achar(9),&
            "x",achar(9),achar(9),achar(9),achar(9),&
            "y",achar(9),achar(9),achar(9),achar(9),&
            "z",achar(9),achar(9),achar(9),achar(9),&
            "v_x",achar(9),achar(9),achar(9),achar(9),&
            "v_y",achar(9),achar(9),achar(9),achar(9),&
            "v_z"

        t = 0.
        do j=1,int(t_max/ H_big)
            t = t + H_big
            ! bodies are uncoupled; faster if integrated separately
            call bulirsch_stoer_step(1, x(i,:), v(i,:), mass(i), H_big)
            
!             print *, 't = '
!             print *, j * H_big
            if (mod(j, int(t_print / H_big)).eq.0) then
                write(save_file_unit,data_format) &
                t, &
                x(i,1), &
                x(i,2), &
                x(i,3), &
                v(i,1), &
                v(i,2), &
                v(i,3)
            endif
        end do
        close(save_file_unit)
    end do


    ! ! ! CHECK FINAL RESULTS
    call state_final%initialize(n_bodies, x, v, mass)

    print *, ""
    print *, 'energy (final): '
    print *,  state_final%energy

    print *, "overall conservation:"
    print *, "energy err:"
    print *, abs( (state_final%energy - state_initial%energy) / state_initial%energy )


end program


