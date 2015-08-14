module functional_mod
    use iso_c_binding, only: dp => c_double
    implicit none
    private
    type, public :: functional
        real(dp) :: energy
        real(dp), allocatable :: grad(:,:,:,:)
        logical :: is_built = .false.
        contains
            procedure :: build
            procedure :: add
    end type
contains
    pure subroutine build(fun, g)
        use grid_mod, only: grid
        type(grid), intent(in) :: g
        class(functional), intent(inout) :: fun
        integer :: err
        if (fun%is_built) then
            return
        else
            allocate( fun%grad(g%no,g%nx,g%ny,g%nz), stat=err, source=0._dp)
            ! if (err /= 0) print *, "fun%grad(g%nx,g%ny,g%nz,g%no): Allocation request denied"
            fun%is_built = .true.
        end if
    end subroutine
    pure subroutine add(fun1, fun2)
        class(functional), intent(inout) :: fun1
        type(functional), intent(in) :: fun2
        fun1%energy = fun1%energy + fun2%energy
        fun1%grad = fun1%grad + fun2%grad
    end subroutine
end module
