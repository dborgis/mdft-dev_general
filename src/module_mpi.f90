module module_mpi

  use precision_kinds, only: dp
  implicit none

  !
  ! MPI header - no idea if it is a good one, the most modern one, etc. That's for testing.
  !
  include 'mpif.h'

  type :: mpi_type
    integer :: nproc, rang, code
  end type
  type(mpi_type) :: mpi

contains

  subroutine init
    implicit none
    call mpi_init(mpi%code)
  end subroutine init

  subroutine finalize
    implicit none
    call mpi_finalize(mpi%code)
  end subroutine finalize

end module module_mpi
