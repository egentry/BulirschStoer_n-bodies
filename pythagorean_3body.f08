program pythagorean_3body
    implicit none
    ! to do:
    interface
        PURE real(kind=16) FUNCTION rational_interp(n_in, x_in, y_in, x_out)
            integer,                         intent(in)  :: n_in
            real(kind=16), dimension(n_in),  intent(in)  :: x_in, y_in
            real(kind=16),                   intent(in)  :: x_out
        END function rational_interp
        
        PURE real(kind=16) FUNCTION calc_energy(n_bodies, x, v, mass)
            integer,                               intent(in)  :: n_bodies
            real(kind=16), dimension(n_bodies,3),  intent(in)  :: x, v
            real(kind=16), dimension(n_bodies),    intent(in)  :: mass
        end function calc_energy

        PURE real(kind=16) FUNCTION calc_angular_momentum_scalar(n_bodies, x, v, mass)
            integer,                               intent(in)  :: n_bodies
            real(kind=16), dimension(n_bodies,3),  intent(in)  :: x, v
            real(kind=16), dimension(n_bodies),    intent(in)  :: mass
        end function calc_angular_momentum_scalar

        subroutine double_array_size_3_index(array_in)
            real(kind=16), allocatable, intent(inout)     :: array_in(:,:,:)
        end subroutine
    end interface

    integer, parameter                       :: n_bodies = 3
    ! 3-position, 3-velocity
    real(kind=16), dimension(n_bodies, 3)    :: x, v
    real(kind=16), dimension(n_bodies)       :: mass

    real(kind=16), parameter     :: H_big = 1./real(10**4, kind=16) ! large step size
    real(kind=16)                :: t
    real(kind=16), parameter     :: t_max = 70.
    real(kind=16)                :: energy_initial, angular_momentum_initial
    real(kind=16)                :: energy,         angular_momentum
    real(kind=16)                :: energy_err,     angular_momentum_err
    real(kind=16), parameter     :: energy_tol           = 1d-11
    real(kind=16), parameter     :: angular_momentum_tol = 1d-6
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

    energy_initial = calc_energy(n_bodies, x, v, mass)
    energy         = energy_initial

    angular_momentum_initial = calc_angular_momentum_scalar(n_bodies, x, v, mass)
    angular_momentum         = angular_momentum_initial

    print *, 'energy (initial): '
    print *,  energy_initial

    print *, 'angular momentum (initial): '
    print *,  angular_momentum_initial

    t = 0.
    do i=1,int(t_max/ H_big)
        t = t + H_big
        call bulirsch_stoer_step(n_bodies, x, v, mass, H_big, energy_tol, angular_momentum_tol)
        
        if (mod(i, int(.1 / H_big)).eq.0) then
!             print *, 't = '
!             print *, i * H_big
            do j=1, n_bodies
                write(save_file_units(j),data_format) t, x(j,1), x(j,2), x(j,3)
            end do
        endif
    end do

    energy           = calc_energy(                 n_bodies, x, v, mass)
    angular_momentum = calc_angular_momentum_scalar(n_bodies, x, v, mass)

    energy_err = abs((energy - energy_initial) / energy_initial)
    angular_momentum_err = abs(angular_momentum - angular_momentum_initial)

    print *, "overall conservation:"
    print *, "energy err:"
    print *, energy_err
    print *, "angular momentum err:"
    print *, angular_momentum_err

end program

subroutine double_array_size_3_index(array_in)
    ! assumes 3 index array
    ! assumes 1st index doubled
    implicit none
    real(kind=16), allocatable, intent(inout)     :: array_in( :,:,:)
    real(kind=16), allocatable                    :: array_tmp(:,:,:)
    integer, dimension(3)                         :: array_in_shape
    integer :: i, j, k  ! loop variables


    array_in_shape = shape(array_in)

    allocate(array_tmp(array_in_shape(1)*2, array_in_shape(2), array_in_shape(3)))

    do concurrent(i=1:array_in_shape(1), j=1:array_in_shape(2), k=1:array_in_shape(3))
        array_tmp(i,j,k) = array_in(i,j,k)
    end do  

    call move_alloc(array_tmp, array_in)

end subroutine double_array_size_3_index

subroutine double_array_size_1_index(array_in)
    ! assumes 1 index array
    ! assumes 1st index doubled
    ! to do:
    !   merge with double_array_size_3_index

    implicit none
    real(kind=16), allocatable, intent(inout)     :: array_in( :)
    real(kind=16), allocatable                    :: array_tmp(:)
    integer, dimension(1)                        :: array_in_shape
    integer :: i  ! loop variables


    array_in_shape = shape(array_in)

    allocate(array_tmp(array_in_shape(1)*2))

    do concurrent(i=1:array_in_shape(1))
        array_tmp(i) = array_in(i)
    end do  

    call move_alloc(array_tmp, array_in)

end subroutine double_array_size_1_index

subroutine bulirsch_stoer_step(n_bodies, x_in, v_in, mass, H_big, energy_tol, angular_momentum_tol)
    ! THINGS TO FIX:
    !   this runs pretty slowly.
    !       -I either need a small H for the entire run,
    !        or a large n_steps_max,
    !        but I think it runs in ~O(n_steps_max^3)
    !        so allowing more n_steps just grinds to a halt

    implicit none
    interface 
        subroutine double_array_size_3_index(array_in)
            ! why do I keep having to declare this interface multiple times?
            real(kind=16), allocatable, intent(inout)     :: array_in( :,:,:)
        end subroutine
        subroutine double_array_size_1_index(array_in)
            ! why do I keep having to declare this interface multiple times?
            real(kind=16), allocatable, intent(inout)     :: array_in(:)
        end subroutine
    end interface

    integer,                               intent(in)     :: n_bodies
    real(kind=16), dimension(n_bodies, 3), intent(inout)  :: x_in, v_in
    real(kind=16), dimension(n_bodies),    intent(in)     :: mass
    real(kind=16),                         intent(in)     :: H_big, energy_tol, angular_momentum_tol

    ! function return types
    real(kind=16) :: rational_interp, calc_energy, calc_angular_momentum_scalar

    real(kind=16)                :: energy_initial, angular_momentum_initial
    real(kind=16)                :: energy,         angular_momentum
    real(kind=16)                :: energy_err,     angular_momentum_err

    ! keep track of results of x,v as a function of n used for modified midpoint
    real(kind=16), allocatable              :: x_tmp(:,:,:), v_tmp(:,:,:)
    real(kind=16), allocatable              :: h_small_array(:)

    ! keep track of highest n used for modified midpoint
    !   and how many 'n' the current x_tmp and v_tmp arrays can fit
    integer                                 :: n_steps
    integer                                 :: n_steps_allocated
    integer, parameter                      :: n_steps_max = 1024

    real(kind=16), parameter  :: h0 = 0. ! for use with richardson extrap.
    !temporary holder of nth modified midpoint result at t = t+H
    real(kind=16), dimension(n_bodies, 3)   :: x_n, v_n

    logical      :: continue_bulirsch_stoer
    integer      :: i,j ! loop variables

    ! set initial values
    continue_bulirsch_stoer  = .True.

    n_steps = 0
    n_steps_allocated = 16

    allocate(x_tmp(n_steps_allocated/2,n_bodies,3))
    allocate(v_tmp(n_steps_allocated/2,n_bodies,3))
    allocate(h_small_array(n_steps_allocated/2))


    energy_initial           = calc_energy(n_bodies, x_in, v_in, mass)
    angular_momentum_initial = calc_angular_momentum_scalar(n_bodies, x_in, v_in, mass)


    do while(continue_bulirsch_stoer.eqv..True.)
        !increase x_tmp, v_tmp arrays if they are full
        if (n_steps .eq. n_steps_allocated) then
            call double_array_size_3_index(x_tmp)
            call double_array_size_3_index(v_tmp)
            call double_array_size_1_index(h_small_array)
            n_steps_allocated = n_steps_allocated * 2
        endif

        n_steps = n_steps + 2

        h_small_array(n_steps/2) = H_big / n_steps

        call modified_midpoint(n_bodies, x_in, v_in, mass, H_big, n_steps, x_n, v_n)
        do concurrent (i=1:n_bodies, j=1:3)
            x_tmp(n_steps/2, i,j) = x_n(i,j)
            v_tmp(n_steps/2, i,j) = v_n(i,j)
        end do

!         energy           = calc_energy(                 n_bodies, x_n, v_n, mass)
!         angular_momentum = calc_angular_momentum_scalar(n_bodies, x_n, v_n, mass)

!         energy_err           = abs((energy - energy_initial) / energy_initial)
!         angular_momentum_err = abs(angular_momentum - angular_momentum_initial)
!         print *, '-- (no_interp)'
!         print *, 'for n = '
!         print *, n_steps
!         print *, 'energy'
!         print *, energy
!         print *, 'energy error: '
!         print *, energy_err
!         print *, 'ang error: '
!         print *, angular_momentum_err
!         print *, '--'

!         print *, 'v_n before: '
!         print *, v_n
!         print *, 'x_n before: '
!         print *, x_n

        do concurrent (i=1:n_bodies, j=1:2)
            x_n(i,j) = rational_interp(n_steps/2, h_small_array(1:n_steps/2), &
                x_tmp(1:n_steps/2, i,j), h0)
            v_n(i,j) = rational_interp(n_steps/2, h_small_array(1:n_steps/2), &
                v_tmp(1:n_steps/2, i,j), h0)
        end do
        do concurrent (i=1:n_bodies)
            ! enforce z=0 plane, otherwise extrapolator can get confused
            x_n(i,3) = 0.
            v_n(i,3) = 0.
        end do

!         print *, 'v_n after: '
!         print *, v_n
!         print *, 'x_n after: '
!         print *, x_n

    
        energy           = calc_energy(                 n_bodies, x_n, v_n, mass)
        angular_momentum = calc_angular_momentum_scalar(n_bodies, x_n, v_n, mass)

        energy_err           = abs((energy - energy_initial) / energy_initial)
        angular_momentum_err = abs(angular_momentum - angular_momentum_initial)

!         print *, '--'
!         print *, 'for n = '
!         print *, n_steps
!         print *, 'energy'
!         print *, energy
!         print *, 'energy error: '
!         print *, energy_err
!         print *, 'ang error: '
!         print *, angular_momentum_err
!         print *, '--'

        if (energy_err.lt.energy_tol) then
            if (angular_momentum_err.lt.angular_momentum_tol) then
                continue_bulirsch_stoer = .False.
!                 print *, 'step ended by finding tolerance'
!                 print *, 'for n = '
!                 print *, n_steps
!                 print *, 'energy err:'
!                 print *, energy_err
!                 print *, 'ang momentum err:'
!                 print *, angular_momentum_err
            end if
        end if
        
        if ( n_steps .ge. n_steps_max ) then
            continue_bulirsch_stoer = .False.
            print *, 'step ended by n_steps cutoff'
            print *, 'for n = '
            print *, n_steps
            print *, 'energy err:'
            print *, energy_err
            print *, 'angular momentum err:'
            print *, angular_momentum_err
        end if

    end do 

    x_in = x_n
    v_in = v_n

end subroutine bulirsch_stoer_step


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

subroutine modified_midpoint(n_bodies, x_in, v_in, mass, H_big, n_steps, x_out, v_out)
    ! follows the algorithm of Eq 16.3.2 in Numerical Recipes, v2
    implicit none
    integer,                               intent(in)    :: n_bodies, n_steps  
    real(kind=16), dimension(n_bodies, 3), intent(in)    :: x_in, v_in
    real(kind=16), dimension(n_bodies   ), intent(in)    :: mass
    real(kind=16),                         intent(in)    :: H_big
    real(kind=16), dimension(n_bodies, 3), intent(out)   :: x_out, v_out

    real(kind=16)  :: h_small   ! fortran is case insensitive; can't just use "H" and "h"

    real(kind=16), dimension(n_bodies, 3) :: force
    real(kind=16), dimension(n_bodies, 3) :: x_previous, v_previous
    real(kind=16), dimension(n_bodies, 3) :: x_current,  v_current
    real(kind=16), dimension(n_bodies, 3) :: x_next,     v_next

    integer :: i, j, m

    h_small = H_big / n_steps

    ! n = 0 base case
    x_previous = x_in   
    v_previous = v_in

    ! n = 1 base case
    call calc_force(n_bodies, x_previous, mass, force)
    do concurrent (i=1:n_bodies, j=1:3)
        x_current(i,j) = x_previous(i,j) + h_small * v_previous(i,j)
        v_current(i,j) = v_previous(i,j) + h_small * force(i,j) / mass(i)
    end do

    do m = 1, n_steps-1
        ! "*_previous" means m-1
        ! "*_current"  means m
        ! "*_next"     means m+1
        call calc_force(n_bodies, x_current, mass, force)
        do concurrent (i=1:n_bodies, j=1:3)
            x_next(i,j) = x_previous(i,j) + 2*h_small * v_current(i,j)
            v_next(i,j) = v_previous(i,j) + 2*h_small * force(i,j) / mass(i)
        end do

        x_previous = x_current
        v_previous = v_current

        x_current = x_next
        v_current = v_next        

    end do

    ! final step
    !   how to match my names to the Eq 16.3.2 subscripts:
    !       z_n     - "current"
    !       z_n-1   - "previous"
    call calc_force(n_bodies, x_current, mass, force)
    do concurrent (i=1:n_bodies, j=1:3)
        x_out(i,j) = 1./2 * ( x_current(i,j) + x_previous(i,j) &
                            + h_small * v_current(i,j) )

        v_out(i,j) = 1./2 * ( v_current(i,j) + v_previous(i,j) &
                            + h_small * force(i,j) / mass(i) )
    end do


end subroutine modified_midpoint

pure function rational_interp(n_in, x_in, y_in, x_out) result(y_out)
    ! rational interpolation, following the alorithm laid out in
    ! Numerical Recipes (v.2), Section 3.2

    ! inputs:
    !   n_in    - integer - number of input data points
    !   x_in    - n_in x 1 float array - locations of observed data
    !   y_in    - n_in x 1 float array - observed data
    !       such that y_i = y(x_i) for x_i in [1, n_in]

    !   x_in    - float - locations of interpolated data

    ! outputs:
    !   y_out -   float - interpolated data
    !       such that y_i = y(x_i) for x_i in [1, n_out]

    ! side effects:
    !   none (enforced)

    ! assumes:
    !   all floats are double precision

    ! features for later:
    !   allow interpolation (/ extrapolation) for a list of values
    !       - key difference will be in how c,d tableaus are handled
    !       - lambda functions would be nice, but slow (?)
    !           and not native in modern Fortran


    implicit none

    integer,                         intent(in)  :: n_in
    real(kind=16), dimension(n_in),  intent(in)  :: x_in, y_in
    real(kind=16),                   intent(in)  :: x_out
    real(kind=16)                                :: y_out

    real(kind=16), dimension(n_in, n_in)         :: c_tableau
    real(kind=16), dimension(n_in, n_in)         :: d_tableau

    ! base case of R
    real(kind=16), parameter                     :: R = 0.
    real(kind=16), dimension(n_in)               :: R_i

    integer     :: m, i     ! loop variable


    ! m = 1 level
    ! NOTE:  m is in [0, n_in-1], but we'll index the m column by m+1
    !   (such that for m=0 base case, it'll be in column 1)
    !   ( will be most confusing for "x_in(i+m)" rather than "x_in(i+m+1)")

    m = 0 
    R_i               = y_in
    c_tableau(:,m+1)  = R_i - R
    d_tableau(:,m+1)  = R_i - R

    do m = 1, n_in-1
        ! Equations 3.2.8 in Numerical Recipes, v2.

        do concurrent (i=1:(n_in-m))
            d_tableau(i, m+1) = c_tableau(i+1,m) &
            * ( c_tableau(i+1,m) - d_tableau(i,m) ) &
            / ( ((x_out - x_in(i)) / (x_out - x_in(i+m)) * d_tableau(i,m)) &
                - c_tableau(i+1,m) )

            c_tableau(i, m+1) = ((x_out - x_in(i)) / (x_out - x_in(i+m))) &
            * d_tableau(i,m) &
            * ( c_tableau(i+1,m) - d_tableau(i,m) ) &
            / ( ((x_out - x_in(i)) / (x_out - x_in(i+m)) * d_tableau(i,m)) &
                - c_tableau(i+1,m) )   
        end do
    end do

    y_out = sum(c_tableau(1,:))

!     !   ! alternatively:
!         y_out = 0.
!         do m=1, n_in
!             y_out = y_out + d_tableau(n_in-(m-1),m)
!         end do
    

end function rational_interp
