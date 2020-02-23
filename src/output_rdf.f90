!===================================================================================================================================
SUBROUTINE output_rdf (array,filename)
!===================================================================================================================================
! Returns the radial distribution functions of each couple solute site - solvent site
! Since the supercell is orthorombic, the maximum range of the RDF is given by half the length of the largest diagonal of a cubic box
! the length of which is the smallest length of our cell.
! The number of bin of the histogram is given by Rice Rule for now. See http://en.wikipedia.org/wiki/Histogram
! A much better (and more complicated) way of chosing this number would be to follow this publication:
! www.proba.jussieu.fr/mathdoc/textes/PMA-721.pdf

! slightly change by Daniel 9-12-2018: print if rdf(bin) /= 0 instead of >0. Works also for negative densities
! (like charge densities)


    use precision_kinds, only: dp
    use module_solute, only: solute
    use module_solvent, only: solvent
    use mathematica, only: deduce_optimal_histogram_properties
    use module_grid, only: grid
    use module_input, only: getinput


    implicit none

    REAL(dp), DIMENSION(grid%nx,grid%ny,grid%nz), INTENT(IN) :: array
    CHARACTER(80), INTENT(IN) :: filename
    REAL(dp) :: RdfMaxRange, dr
    REAL(dp), ALLOCATABLE :: rdf(:), solute_site_gr(:,:)
    INTEGER:: n,bin,nbins
    TYPE :: errortype
        LOGICAL :: found
        CHARACTER(180) :: msg
    END TYPE
    TYPE (errortype) :: error
    real(dp) :: xbin
    logical :: dontPrintZeros
    character(len=6) :: string, string2
    character(80) :: filename_xvg, filename_out

    rdfmaxrange = minval(grid%length)/2._dp
    ! we dont use this anymore. Not suited to our "powder averaging" kind of grid results.
    ! CALL deduce_optimal_histogram_properties( product(grid%n_nodes), rdfmaxrange, nbins, dr )
    dr = getinput%dp( "rdf_dr", defaultvalue=0.1_dp, assert=">0") !Angstrom
    nbins = int( rdfmaxrange/dr ) +1

    allocate (rdf(nbins), source=0._dp)
    allocate( solute_site_gr(nbins, size(solute%site)), source=0.0_dp )

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
     filename_xvg = trim(adjustl(filename))//'.xvg'
     filename_out = trim(adjustl(filename))//'.out'

    open (10, file=filename_xvg) ! filename is intent(in), typically "output/rdf.out"
    open(11, file=filename_out)
    write(*,*)
    write(10,*) '@ xaxis label "Distance (\cE\C)"'
    write(10,*) '@ yaxis label "Radial distribution function"'
    write(10,*) '@ TYPE xy'
    write(10,*) '@ legend on'
    write(10,*) '@ legend box on'
    write(10,*) '@ legend loctype view'
    write(10,*) '@ legend 0.78, 0.8'
    write(10,*) '@ legend length 4'
    do n=1, size(solute%site) ! loop over all sites of the solute
        call histogram_3d (array(:,:,:), solute%site(n)%r, rdf)
        write(string,'(i5)') n-1
        write(string2,'(i5)') n
        write(10,*) '@ s'//trim(adjustl(string))//' legend "solute site '//trim(adjustl(string2))//'"'
        dontPrintZeros = .false.
        do bin=1,nbins
            xbin = real((bin-0.5)*dr) ! we use the coordinates of the middle of the bin. first bin from x=0 to x=dr is written has dr/2. 2nd bin [dr,2dr] has coordinate 1.5dr
            solute_site_gr(bin,n) = rdf(bin)
            if( .not. dontPrintZeros ) then ! print the zeros
                write(10,*) xbin, rdf(bin) ! For bin that covers 0<r<dr, I print at 0.5dr, i.e., at the middle of the bin
                if( rdf(bin)/=0 ) dontPrintZeros = .true.
            else if( dontPrintZeros ) then
                if( rdf(bin)/=0 ) then
                    write(10,*) xbin, rdf(bin)
                end if
            end if
        end do
    end do
    close(10)

   do bin=1, nbins
      xbin = real((bin-0.5)*dr)
      write(11,'(11F7.3)') xbin,(solute_site_gr(bin,n), n=1,min(10,size(solute%site)))  !limited to 10 sites
   end do
   close(11)

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
  integer :: ix, iy
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

            deallocate( occurrence )

        end subroutine histogram_3d

END SUBROUTINE output_rdf
