module density_mod
    use iso_c_binding, only: dp => c_double
    implicit none
    private
    type, public :: density
        real(dp), allocatable :: rho(:,:,:,:)
        real(dp) :: bulk = 0.0333_dp
        logical :: is_built = .false.
        contains
        procedure :: build
        procedure, nopass :: reshapeto1d
        procedure, nopass :: reshapeto4d
        procedure :: is_negative_somewhere
    end type
contains
    subroutine build(den, gr)
        use input, only: getinput%dp
        use grid_mod, only: grid
        implicit none
        class(density), intent(inout) :: den
        type(grid) :: gr
        integer :: nx, ny, nz, no
        if (den%is_built) then
            return
        else
            nx = grid%nx
            ny = grid%ny
            nz = grid%nz
            no = grid%no
            den%bulk = getinput%dp("bulk_density", defaultvalue=0.0333_dp)
            allocate( den%rho(no,nx,ny,nz) ,source=den%bulk)
            den%is_built = .true.
        end if
    end subroutine
    subroutine reshapeto1d
    end subroutine
    subroutine reshapeto4d
    end subroutine
    function is_negative_somewhere(den)
        class(density), intent(in) :: den
        logical :: is_negative_somewhere
        if (any(den%rho<0._dp)) then
            is_negative_somewhere = .true.
        else
            is_negative_somewhere = .false.
        end if
    end function
end module
