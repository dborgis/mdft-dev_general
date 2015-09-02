module external_potential_mod
    use iso_c_binding, only: dp => c_double
    implicit none
    private
    type, public :: external_potential
        real(dp), allocatable :: tot(:,:,:,:)
        logical :: is_built = .false.
        contains
        procedure :: build
    end type
contains
    subroutine build(vxt, gr)
        use grid_mod, only: grid
        class(external_potential), intent(inout) :: vxt
        type(grid), intent(in) :: gr
        integer :: err
        if (vxt%is_built) return
        if (.not. grid%is_built) error stop "gr is not built in external_potential_mod"
        allocate( vxt%tot(grid%no, grid%nx, grid%ny, grid%nz), stat=err, source=0._dp)
        if (err /= 0) print *, "vxt: Allocation request denied"
        call build_lj
        vxt%is_built = .true.
    end subroutine
    subroutine build_lj
    end subroutine
end module
