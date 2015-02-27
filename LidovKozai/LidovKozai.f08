program LidovKozai
    use physics
    use bulirsch_stoer
    implicit none
    ! to do:
    !   figure out how to make it so I don't have to explicitly use 'accuracy.f08'
    !       within bulirsch_stoer.f08
    !       -- that would help make this a bit easier to use modularly

    ! the main things you need to set are:
    !    In physics_LidovKozai.f08:
    !       - precision - sets float precision (e.g. 8 = double)
    !       - energy_tol - fractional tolerance for energy conservation *each step*
    !       - Body masses
    !           - M_star    - innermost star
    !           - M_planet  - inner companion to star; should be low mass
    !           - M_outer   - outer, long-period binary companion
    !       - Semi-major axes (unperturbed)
    !           - a_planet_0 - Planet - Star binary
    !           - a_outer_0  - outer object - Star binary
    !       - Theta
    !           - conserved integral of motion
    !           - Theta = (1 - e_planet**2) cos i
    !       - x_0 (= 1- e_0^2)
    !           - sets initial eccentricity, 
    !             which determines inclination, argument of periastron
    !

    !   In this file:
    !   
    !       - t_max
    !           - how long should the integration continue?
    !       - t_print
    !           - how often should you print?

 


    integer, parameter                              :: n_bodies = 3
    ! 3-position, 3-velocity
    real(kind=precision), dimension(n_bodies, 3)    :: x, v
    real(kind=precision), dimension(n_bodies)       :: mass
    real(kind=precision)                            :: mass_total

    real(kind=precision), parameter     :: H_big = 1./real(10**(3), kind=precision) ! large step size
    real(kind=precision)                :: t
    real(kind=precision), parameter     :: t_max   = 100000.
    real(kind=precision), parameter     :: t_print = .1  ! print every t_print,

    ! for transforming into Center of Mass Frame
    real(kind=precision)                :: v_x_COM, v_y_COM, v_z_COM ! to remove center of mass velocity
    real(kind=precision)                :: x_COM, y_COM, z_COM ! to remove center of mass position


    ! defined in the physics module -- changes for each potential
    type(accuracy_state)                :: state_initial, state_final

    integer                             :: i,j,k        ! loop variables

    integer, dimension(n_bodies)        :: save_file_units
    character(len=6)                    :: jth_body_str ! tmp variable used in filenames
    character(len=15)                   :: filename     ! tmp variable used in filenames
    character(len=15)                   :: data_format  ! for saving

    do j=1, n_bodies
        write (jth_body_str, "(I3.3)") j
        filename = 'body_'//trim(jth_body_str)//'.dat'

        open(newunit=save_file_units(j), file=filename, status='unknown')
        write(save_file_units(j), *) &
            "t",achar(9),achar(9),achar(9),achar(9),&
            "x",achar(9),achar(9),achar(9),achar(9),&
            "y",achar(9),achar(9),achar(9),achar(9),&
            "z",achar(9),achar(9),achar(9),achar(9),&
            "v_x",achar(9),achar(9),achar(9),achar(9),&
            "v_y",achar(9),achar(9),achar(9),achar(9),&
            "v_z"
    end do
    data_format = "(7e13.5)"

    ! INITIALIZE

    mass(1) = M_star
    mass(2) = M_planet
    mass(3) = M_outer

    mass(1) = M_star
    x(1, 1) =  0.                        ! body 1, x position   
    x(1, 2) =  0.                        ! body 1, y position
    x(1, 3) =  0.                        ! body 1, z position
    v(1, 1) =  0.                        ! body 1, x velocity 
    v(1, 2) =  0.                        ! body 1, y velocity 
    v(1, 3) =  0.                        ! body 1, z velocity


    mass(2) =  M_planet
    x(2, 1) =  r_planet_0 * sin(i_0)     ! body 2, x position 
    x(2, 2) =  0.                        ! body 2, y position 
    x(2, 3) =  r_planet_0 * cos(i_0)     ! body 2, z position 
    v(2, 1) =  0.                        ! body 2, x velocity 
    v(2, 2) =  v_planet_0                ! body 2, y velocity 
    v(2, 3) =  0.                        ! body 2, z velocity


    mass(3) =  M_outer
    x(3, 1) =  r_outer_0                 ! body 3, x position 
    x(3, 2) =  0.                        ! body 3, y position 
    x(3, 3) =  0.                        ! body 3, z position
    v(3, 1) =  0.                        ! body 3, x velocity 
    v(3, 2) =  v_outer_0                 ! body 3, y velocity 
    v(3, 3) =  0.                        ! body 3, z velocity

    mass_total = sum(mass)
    !! Go into Center of Mass Frame
    v_x_COM = 0.
    v_y_COM = 0.
    v_z_COM = 0.
    x_COM   = 0.
    y_COM   = 0.
    z_COM   = 0.

    do concurrent (i=1:n_bodies)
        v_x_COM = v_x_COM + v(i,1)*mass(i)/mass_total
        v_y_COM = v_y_COM + v(i,2)*mass(i)/mass_total
        v_z_COM = v_z_COM + v(i,3)*mass(i)/mass_total
          x_COM =   x_COM + x(i,1)*mass(i)/mass_total
          y_COM =   y_COM + x(i,2)*mass(i)/mass_total
          z_COM =   z_COM + x(i,3)*mass(i)/mass_total
    end do
    do concurrent (i=1:n_bodies)
        v(i,1) = v(i,1) - v_x_COM
        v(i,2) = v(i,2) - v_y_COM
        v(i,3) = v(i,3) - v_z_COM

        x(i,1) = x(i,1) -   x_COM
        x(i,2) = x(i,2) -   y_COM
        x(i,3) = x(i,3) -   z_COM
    end do


    call state_initial%initialize(n_bodies, x, v, mass)

    print *, 'energy (initial): '
    print *,  state_initial%energy

!     print *, "a_planet_0"
!     print *, a_planet_0
!     print *, r_planet_0
!     print *, e_0 
!     print *, x
!     print *, v

    ! ! ! BEGIN INTEGRATIONS
    t = 0.
    do i=1,int(t_max/ H_big)
        t = t + H_big
        call bulirsch_stoer_step(n_bodies, x, v, mass, H_big)
        
        if (mod(i, int(t_print / H_big)).eq.0) then
!             print *, 't = '
!             print *, i * H_big
            do j=1, n_bodies
                write(save_file_units(j),data_format) &
                t, &
                x(j,1), &
                x(j,2), &
                x(j,3), &
                v(j,1), &
                v(j,2), &
                v(j,3)
            end do
        endif
    end do

    close(save_file_units(1))
    close(save_file_units(2))
    close(save_file_units(3))

    ! ! ! CHECK FINAL RESULTS
    call state_final%initialize(n_bodies, x, v, mass)

    print *, ""
    print *, 'energy (final): '
    print *,  state_final%energy

    print *, "overall conservation:"
    print *, "energy err:"
    print *, abs( (state_final%energy - state_initial%energy) / state_initial%energy )


end program


