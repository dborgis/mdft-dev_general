!===================================================================================================================================
SUBROUTINE output_rdf (array,filename)
!===================================================================================================================================
! Returns the radial distribution functions of each couple solute site - solvent site
! Since the supercell is orthorombic, the maximum range of the RDF is given by half the length of the largest diagonal of a cubic box
! the length of which is the smallest length of our cell.
! The number of bin of the histogram is given by Rice Rule for now. See http://en.wikipedia.org/wiki/Histogram
! A much better (and more complicated) way of chosing this number would be to follow this publication:
! www.proba.jussieu.fr/mathdoc/textes/PMA-721.pdf

    use precision_kinds, only: dp
    use module_solute, only: solute
    use module_solvent, only: solvent
    use mathematica, only: deduce_optimal_histogram_properties
    use module_grid, only: grid

    implicit none

    REAL(dp), DIMENSION(grid%nx,grid%ny,grid%nz), INTENT(IN) :: array
    CHARACTER(50), INTENT(IN) :: filename
    REAL(dp) :: RdfMaxRange, dr
    REAL(dp), ALLOCATABLE :: rdf(:)
    INTEGER:: n,bin,nbins
    TYPE :: errortype
        LOGICAL :: found
        CHARACTER(180) :: msg
    END TYPE
    TYPE (errortype) :: error
    real(dp) :: xbin
    logical :: dontPrintZeros

    if (solvent(1)%nspec/=1) error stop "compute_rdf.f90 is written for 1 solvent species only."

    rdfmaxrange = minval(grid%length)/2._dp
    ! we dont use this anymore. Not suited to our "powder averaging" kind of grid results.
    CALL deduce_optimal_histogram_properties( product(grid%n_nodes), rdfmaxrange, nbins, dr )
    dr = 0.01 ! Angstrom
    nbins = int( rdfmaxrange/dr ) +1

    allocate (rdf(nbins), source=0._dp)
    call UTest_histogram_3D

    !
    ! Compute and print the site-site radial distribution of each solute site
    ! to the grid nodes. Grid nodes are for instance oxygen sites of SPCE water
    ! since in input/solvent.in, oxygen sites are located at {0,0,0}.
    !
    ! About: dontPrintZeros.
    ! We have empty bins in the RDF:
    ! - the ones in the core of the solute, for which we want to print rdf(x in core)=0.
    ! - the ones that are empty because of lack of information, for x > beginning of the first peak.
    ! These last ones we dont want to print.
    ! We thus first print the zero, but as soon as we detect the rdf starts to be nonzero,
    ! we dont print the zeros anymore.
    !
    !
    open (10, file=filename) ! filename is intent(in), typically "output/rdf.out"
    do n=1, size(solute%site) ! loop over all sites of the solute
        call histogram_3d (array(:,:,:), solute%site(n)%r, rdf)
        write(10,*)'# solute site', n
        dontPrintZeros = .false.
        do bin=1,nbins
            xbin = real((bin-0.5)*dr) ! we use the coordinates of the middle of the bin. first bin from x=0 to x=dr is written has dr/2. 2nd bin [dr,2dr] has coordinate 1.5dr
            if( .not. dontPrintZeros ) then ! print the zeros
                write(10,*) xbin, rdf(bin) ! For bin that covers 0<r<dr, I print at 0.5dr, i.e., at the middle of the bin
                if( rdf(bin)>0 ) dontPrintZeros = .true.
            else if( dontPrintZeros ) then
                if( rdf(bin)>0 ) write(10,*) xbin, rdf(bin)
            end if
        end do
        write(10,*) ! skip a line between each solute site so that it appears nicely in xmgrace.
    end do
    close(10)

contains

pure subroutine histogram (Y, X, H)
  ! Y is an array of size N, it contains N values of type real one wants to build an histogram from.
  ! X contains the bounds of the histogram. X(1) is the lower bound.
  ! H contains the histogram. if  x(i)<=y(j)<(i+1)  then h(i) is incremented by 1.
  implicit none
  real, intent(in) :: Y(:)
  real, intent(in) :: X(:)
  integer, intent(out) :: H(:)
  integer :: sizeY, sizeX, sizeH ! size of arrays Y, X and H
  integer :: ix, iy, iz
  sizeY = size(Y)
  sizeX = size(X)
  sizeH = size(H)
  H = 0
  Yvalues: do iY = 1, sizeY
    do iX = 1, sizeX-1
      if( X(iX) <= Y(iY) .and. Y(iY) < X(iX+1) ) then
        H(iX) = H(iX)+1
        cycle Yvalues
      end if
    end do
  end do Yvalues
end subroutine histogram

        !===========================================================================================================================
        SUBROUTINE histogram_3d (data, origin, rdf)
        !===========================================================================================================================
            implicit none
            REAL(dp), intent(in) :: data(grid%nx,grid%ny,grid%nz)
            REAL(dp), intent(out) :: rdf(nbins)
            real(dp), intent(in) :: origin(3)
            INTEGER, allocatable :: occurrence(:)
            REAL(dp) :: r,rsq,xisq,yisq,zisq
            INTEGER :: ibin, ix, iy, iz

            allocate( occurrence (nbins) ,source=0)
            rdf = 0._dp

            ! Transform array(position) in rdf(radialdistance)
            ! counts the total number of appearence of a value in each bin
            do ix=1,grid%nx
                xisq=((ix-1)*grid%dx-origin(1))**2
                do iy=1,grid%ny
                    yisq=((iy-1)*grid%dy-origin(2))**2
                    do iz=1,grid%nz
                        zisq=((iz-1)*grid%dz-origin(3))**2
                        rsq=xisq+yisq+zisq
                        r=sqrt(rsq)
                        ibin=int(r/dr) +1
                        if (ibin>nbins) then
                            cycle
                        else if (ibin<=0) then
                            error stop "I have a negative or null index for the bin in output_rdf.f90"
                        else
                            occurrence(ibin)=occurrence(ibin) +1
                            rdf(ibin) =rdf(ibin) +data(ix,iy,iz)
                        end if
                    end do
                end do
            end do

            if (all(occurrence==0)) then
                print*, "In output_rdf, the histogram is completely empty. all bins have zero occurence"
                error stop
            end if

            ! normalize the rdf
            where (occurrence/=0)
                rdf = rdf/REAL(occurrence,dp)
            elsewhere
                rdf = 0.0_dp
            end where

        end subroutine histogram_3d

        SUBROUTINE UTest_histogram_3D
            IMPLICIT NONE
            REAL(dp), ALLOCATABLE :: nullarray3D (:,:,:)
            REAL(dp), ALLOCATABLE :: rdfUT(:)
            real(dp), parameter :: zerodp=0._dp
            real(dp), parameter :: zerodp3(3)=[zerodp,zerodp,zerodp]
            real(dp), parameter :: epsdp=epsilon(1._dp)

            ALLOCATE (nullarray3D (grid%nx,grid%ny,grid%nz) ,SOURCE=0._dp)
            ALLOCATE (rdfUT(nbins),source=0._dp)

            ! Test 1
            nullarray3D = 0._dp
            CALL histogram_3d (nullarray3D, zerodp3, rdfUT)
            IF (ANY(abs(rdfUT)>epsdp)) error stop "Test 1 in UTest_histogram_3D not passed."

            ! Test 2
            nullarray3D = 100._dp
            CALL histogram_3d (nullarray3D, zerodp3, rdfUT)
            IF (ANY(rdfUT<0._dp)) STOP "Test 2 in UTest_histogram_3D not passed."

            ! Test 3
            nullarray3D = -1._dp
            CALL histogram_3d (nullarray3D, zerodp3, rdfUT)
            IF (ANY(rdfUT>0._dp)) STOP "Test 3 in UTest_histogram_3D not passed."

            ! Test 4
            nullarray3D = 100._dp
            CALL histogram_3d (nullarray3D, zerodp3, rdfUT)
            IF (ANY(rdfUT>100._dp)) STOP "Test 4 in UTest_histogram_3D not passed."

            DEALLOCATE (nullarray3D, rdfUT)
        END SUBROUTINE UTest_histogram_3D

END SUBROUTINE output_rdf
