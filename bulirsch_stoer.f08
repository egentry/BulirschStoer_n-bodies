  module bulirsch_stoer 
    use physics

contains

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

    subroutine bulirsch_stoer_step(n_bodies, x_in, v_in, mass, H_big)
        ! THINGS TO FIX:
        !   this runs pretty slowly.
        !       -I either need a small H for the entire run,
        !        or a large n_steps_max,
        !        but I think it runs in ~O(n_steps_max^3)
        !        so allowing more n_steps just grinds to a halt

        implicit none
        integer,                               intent(in)     :: n_bodies
        real(kind=16), dimension(n_bodies, 3), intent(inout)  :: x_in, v_in
        real(kind=16), dimension(n_bodies),    intent(in)     :: mass
        real(kind=16),                         intent(in)     :: H_big

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

        type(accuracy_state)        :: state_initial, state_current
        call state_initial%initialize(n_bodies, x_in, v_in, mass)


        ! set initial values
        continue_bulirsch_stoer  = .True.

        n_steps = 0
        n_steps_allocated = 16

        allocate(x_tmp(n_steps_allocated/2,n_bodies,3))
        allocate(v_tmp(n_steps_allocated/2,n_bodies,3))
        allocate(h_small_array(n_steps_allocated/2))



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

            call state_current%initialize(n_bodies, x_n, v_n, mass)

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

            if (accuracy_check(state_initial, state_current).eqv..True.) then
                continue_bulirsch_stoer = .False.
!                 print *, 'step ended by finding tolerance'
!                 print *, 'for n = '
!                 print *, n_steps
!                 print *, 'energy err:'
!                 print *, energy_err
!                 print *, 'ang momentum err:'
!                 print *, angular_momentum_err
            end if
            
            if ( n_steps .ge. n_steps_max ) then
                continue_bulirsch_stoer = .False.
                print *, 'step ended by n_steps cutoff'
                print *, 'for n = '
                print *, n_steps
!                 print *, 'energy err:'
!                 print *, energy_err
!                 print *, 'angular momentum err:'
!                 print *, angular_momentum_err
            end if

        end do 

        x_in = x_n
        v_in = v_n

    end subroutine bulirsch_stoer_step

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

end module