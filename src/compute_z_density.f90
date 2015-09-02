!===================================================================================================================================
SUBROUTINE compute_z_density (array, filename)
!===================================================================================================================================
! Compute average density in the plane z {001}, that is the one perpendicular to plane (xy)
  use precision_kinds ,only: dp, i2b
  use system          ,only: spacegrid
  use mathematica     ,only: chop
  use input           ,only: verbose
  implicit none
  type xy
    real(dp) :: x,y
  end type
  type(xy)                  :: nz(spacegrid%n_nodes(3))
  real(dp)                  :: dx
  real(dp), intent(in)      :: array(spacegrid%n_nodes(1),spacegrid%n_nodes(2),spacegrid%n_nodes(3))
  integer(i2b)              :: k,kmax,nodesinxy
  character(50), intent(in) :: filename
  dx=spacegrid%dl(3)
  kmax=spacegrid%n_nodes(3)
  nodesinxy=product(spacegrid%n_nodes(1:2)) ! product of the number of nodes in x and y directions
  ! Compute mean density over x and y
  do k=1,kmax
    nz(k)%x=(k-1)*dx
    nz(k)%y=chop( sum(array(:,:,k))/nodesinxy )
  end do

  open(10, file=filename)
  do k=1,kmax
    write(10,*) nz(k)%x, nz(k)%y
  end do
  close(10)
  if (verbose) print*,"Written ", trim(adjustl(filename))

end subroutine compute_z_density
!===================================================================================================================================
