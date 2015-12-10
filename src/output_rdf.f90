!===================================================================================================================================
SUBROUTINE output_rdf (array,filename)
!===================================================================================================================================
! Returns the radial distribution functions of each couple solute site - solvent site
! Since the supercell is orthorombic, the maximum range of the RDF is given by half the length of the largest diagonal of a cubic box
! the length of which is the smallest length of our cell.
! The number of bin of the histogram is given by Rice Rule for now. See http://en.wikipedia.org/wiki/Histogram
! A much better (and more complicated) way of chosing this number would be to follow this publication:
! www.proba.jussieu.fr/mathdoc/textes/PMA-721.pdf

    use precision_kinds, only: dp, i2b, sp
    use module_solute, only: solute
    use module_solvent, only: solvent
    use mathematica, only: deduce_optimal_histogram_properties, chop
    use module_grid, only: grid

    implicit none

    REAL(dp), DIMENSION(grid%nx,grid%ny,grid%nz), INTENT(IN) :: array
    CHARACTER(50), INTENT(IN) :: filename
    REAL(dp) :: RdfMaxRange, dr
    REAL(dp), ALLOCATABLE :: rdf(:)
    INTEGER(i2b):: n,bin,nbins
    TYPE :: errortype
        LOGICAL :: found
        CHARACTER(180) :: msg
    END TYPE
    TYPE (errortype) :: error
    real(sp) :: xbin
    ! integer, parameter :: sxs = 1000
    ! real(dp) :: xs(sxs) , rdfs(sxs), lastx

    if (.not. allocated(solvent)) then
        print*, "In output_rdf.f90, solvent derived type is not allocated"
        error stop
    end if
    if (solvent(1)%nspec/=1) error stop "compute_rdf.f90 is written for 1 solvent species only."

    rdfmaxrange = minval(grid%length)/2._dp
    CALL deduce_optimal_histogram_properties( product(grid%n_nodes), rdfmaxrange, nbins, dr )

    allocate (rdf(nbins), source=0._dp)
    call UTest_histogram_3D
    !
    ! lastx = (real(nbins,dp)-0.5_dp)*dr
    ! do i = 1, sxs
    !     xs(i) = real(i-1)/sxs * lastx
    ! end do

    open (10, file=filename)
    ! open (12, file=trim(adjustl(filename))//"-spline")
    !
    ! Compute and print histogram (the rdf) for each solute site
    !
    do n=1, size(solute%site) ! loop over all sites of the solute
        call histogram_3d (array(:,:,:), solute%site(n)%r, rdf)

        ! write to output/rdf.out
        write(10,*)'# solute site', n
        write(10,*) 0., 0. ! we impose
        do bin=1,nbins
            xbin = real((bin-0.5)*dr) ! we use the coordinates of the middle of the bin. first bin from x=0 to x=dr is written has dr/2. 2nd bin [dr,2dr] has coordinate 1.5dr
            write(10,*) xbin, real(chop(rdf(bin)))! For bin that covers 0<r<dr, I print at 0.5dr, i.e., at the middle of the bin
        end do
        write(10,*)

        ! call akima_spline( nbins, x, rdf, size(xs), xs, rdfs )
        !
        ! write(12,*) "#solute site", n
        ! do i = 1, size(xs)
        !     write(12,*) real([ xs(i), chop(rdfs(i))  ])
        ! end do
        ! write(12,*)
    end do
    close(10)
    ! close(12)

contains


        !===========================================================================================================================
        SUBROUTINE histogram_3d (data, origin, rdf)
        !===========================================================================================================================
            implicit none
            REAL(dp), intent(in) :: data(grid%nx,grid%ny,grid%nz)
            REAL(dp), intent(out) :: rdf(nbins)
            real(dp), intent(in) :: origin(3)
            INTEGER(i2b), allocatable :: occurrence(:)
            REAL(dp) :: r,rsq,xisq,yisq,zisq
            INTEGER(i2b) :: ibin, ix, iy, iz

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
