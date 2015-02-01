program non_axisymmetric_potential
    use physics
    use bulirsch_stoer
    implicit none
    ! to do:
    !   figure out how to make it so I don't have to explicitly use 'accuracy.f08'
    !       within bulirsch_stoer.f08
    !       -- that would help make this a bit easier to use modularly

    ! the main things you need to set are:
    !    In physics_non_axisymmetric.f08:
    !       - precision - sets float precision (e.g. 8 = double)
    !       - energy_tol - fractional tolerance for energy conservation *each step*
    !       - eccentricity - eccentricity of the potential well
    !       - energy_initial - defines family of orbits
    !

    !   In this file:
    !       - n_bodies
    !           - will automatically start n_bodies of integration,
    !             spaced throughout the initially allowed x_values ( 0 < x_0 < x_max )
    !       - t_max
    !           - how long should the integration continue?
    !       - t_print
    !           - how often should you print?




    integer, parameter                              :: n_bodies = 100
    ! 3-position, 3-velocity
    real(kind=precision), dimension(n_bodies, 3)    :: x, v
    real(kind=precision), dimension(n_bodies)       :: mass

    real(kind=precision), parameter     :: H_big = 1./real(10**(4), kind=precision) ! large step size
    real(kind=precision)                :: t
    real(kind=precision), parameter     :: t_max   = 100.
    real(kind=precision), parameter     :: t_print = .001. ! print every t_print,


    type(accuracy_state)                :: state_initial, state_final

    integer                             :: i,j,k        ! loop variables

    integer, dimension(n_bodies)        :: save_file_units
    character(len=6)                    :: jth_body_str ! tmp variable used in filenames
    character(len=15)                   :: filename     ! tmp variable used in filenames
    character(len=15)                   :: data_format  ! for saving

    ! INITIALIZE

    ! initialize x_0 to a range of 0 < x_0 < x_max
    do concurrent (i=1:n_bodies)
        x(i,1) = x_max * i / (n_bodies + 1)
        x(i,2) = 0
        x(i,3) = 0
    end do

    mass = 1.

    call calc_initial_velocity(n_bodies, x, v)


    do j=1, n_bodies
        write (jth_body_str, "(I3.3)") j
        filename = 'body_'//trim(jth_body_str)//'.dat'

        open(newunit=save_file_units(j), file=filename, status='unknown')
        write(save_file_units(j), *) "t",achar(9),achar(9),achar(9),achar(9),&
            "x",achar(9),achar(9),achar(9),achar(9),&
            "y",achar(9),achar(9),achar(9),achar(9),&
            "z",achar(9),achar(9),achar(9),achar(9),&
            "v_x",achar(9),achar(9),achar(9),achar(9),&
            "v_y",achar(9),achar(9),achar(9),achar(9),&
            "v_z"
    end do
    data_format = "(7e13.5)"


    call state_initial%initialize(n_bodies, x, v, mass)

    print *, 'energy (initial): '
    print *,  state_initial%energy

    print *, '---'
    print *, 'x_max:'
    print *, x_max
    print *, 'x_0: '
    print *, x(:,1)

    print *, ''
    print *, 'v_y: '
    print *, v(:,2)


    ! ! ! BEGIN INTEGRATIONS
    do i = 1,n_bodies
        t = 0.
        do j=1,int(t_max/ H_big)
            t = t + H_big
            ! bodies are uncoupled; faster if integrated separately
            call bulirsch_stoer_step(1, x(i,:), v(i,:), mass(i), H_big)
            
!             print *, 't = '
!             print *, j * H_big
            if (mod(j, int(t_print / H_big)).eq.0) then
                write(save_file_units(i),data_format) t, x(i,1), x(i,2), x(i,3), v(i,1), v(i,2), v(i,3)
            endif
        end do
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


