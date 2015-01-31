program pythagorean_3body
    use physics
    use bulirsch_stoer
    implicit none
    ! to do:
    !   figure out how to make it so I don't have to explicitly use 'accuracy.f08'
    !       within bulirsch_stoer.f08
    !       -- that would help make this a bit easier to use modularly

    integer, parameter                       :: n_bodies = 3
    ! 3-position, 3-velocity
    real(kind=16), dimension(n_bodies, 3)    :: x, v
    real(kind=16), dimension(n_bodies)       :: mass

    real(kind=16), parameter     :: H_big = 1./real(10**(4), kind=16) ! large step size
    real(kind=16)                :: t
    real(kind=16), parameter     :: t_max = 70.

    type(accuracy_state)         :: state_initial, state_final

    integer                      :: i,j     ! loop variables

    integer, dimension(n_bodies) :: save_file_units
    character(len=6)             :: jth_body_str ! tmp variable used in filenames
    character(len=15)            :: filename     ! tmp variable used in filenames
    character(len=15)            :: data_format  ! for saving

    do j=1, n_bodies
        write (jth_body_str, "(I3.3)") j
        filename = 'body_'//trim(jth_body_str)//'.dat'

        open(newunit=save_file_units(j), file=filename, status='unknown')
        write(save_file_units(j), *) "t",achar(9),achar(9),achar(9),achar(9),&
            "x",achar(9),achar(9),achar(9),achar(9),&
            "y",achar(9),achar(9),achar(9),achar(9),&
            "z"
    end do
    data_format = "(4e13.5)"


    !initialize
    mass(1) = 3.
    x(1, 1) = 1.          ! body 1, x position   
    x(1, 2) = 3.          ! body 1, y position
    x(1, 3) = 0.          ! body 1, z position


    mass(2) =  4.
    x(2, 1) = -2.          ! body 2, x position 
    x(2, 2) = -1.          ! body 2, y position 
    x(2, 3) =  0.          ! body 2, z position 

    mass(3) =  5.
    x(3, 1) =  1.          ! body 3, x position 
    x(3, 2) = -1.          ! body 3, y position 
    x(3, 3) =  0.          ! body 3, z position

    v = 0.     ! all bodies, all velocities

    call state_initial%initialize(n_bodies, x, v, mass)

    print *, 'energy (initial): '
    print *,  state_initial%energy

    print *, 'angular momentum (initial): '
    print *,  state_initial%angular_momentum

    t = 0.
    do i=1,int(t_max/ H_big)
        t = t + H_big
        call bulirsch_stoer_step(n_bodies, x, v, mass, H_big)
        
        if (mod(i, int(.1 / H_big)).eq.0) then
!             print *, 't = '
!             print *, i * H_big
            do j=1, n_bodies
                write(save_file_units(j),data_format) t, x(j,1), x(j,2), x(j,3)
            end do
        endif
    end do

    call state_final%initialize(n_bodies, x, v, mass)

    print *, ""
    print *, 'energy (final): '
    print *,  state_final%energy

    print *, 'angular momentum (final): '
    print *,  state_final%angular_momentum

    print *, "overall conservation:"
    print *, "energy err:"
    print *, abs( (state_final%energy - state_initial%energy) / state_initial%energy )
    print *, "angular momentum err:"
    print *, abs( state_final%angular_momentum - state_initial%angular_momentum )

end program


