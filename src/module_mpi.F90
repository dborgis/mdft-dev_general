module module_mpi

  !
  ! MPI header - no idea if it is a good one, the most modern one, etc. That's for testing.
  !
#ifndef MPI_FORTRAN_HAS_MODULE
  include 'mpif.h'
#else
  use mpi
#endif
  
  implicit none

  type :: impi_type
    integer :: nproc, rang, code
  end type impi_type
  type(impi_type) :: impi

contains

  subroutine init
    implicit none
    call mpi_init(impi%code)
  end subroutine init

  subroutine finalize
    implicit none
    call mpi_finalize(impi%code)
  end subroutine finalize

end module module_mpi
