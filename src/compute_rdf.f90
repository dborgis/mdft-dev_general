!===================================================================================================================================
SUBROUTINE compute_rdf (array,filename)
!===================================================================================================================================
! Returns the radial distribution functions of each couple solute site - solvent site
! Since the supercell is orthorombic, the maximum range of the RDF is given by half the length of the largest diagonal of a cubic box
! the length of which is the smallest length of our cell.
! The number of bin of the histogram is given by Rice Rule for now. See http://en.wikipedia.org/wiki/Histogram
! A much better (and more complicated) way of chosing this number would be to follow this publication:
! www.proba.jussieu.fr/mathdoc/textes/PMA-721.pdf

    USE precision_kinds ,ONLY: dp, i2b, sp
    USE system          ,ONLY: nb_species, spaceGrid, soluteSite, solventSite
    
    IMPLICIT NONE

    REAL(dp), DIMENSION(spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3)), INTENT(IN) :: array
    CHARACTER(50), INTENT(IN) :: filename
    REAL(dp) :: RdfMaxRange, dr, r
    REAL(dp), ALLOCATABLE :: rdf(:)
    INTEGER(i2b):: n,bin,nbins
    TYPE :: errortype
        LOGICAL :: found
        CHARACTER(180) :: msg
    END TYPE
    TYPE (errortype) :: error
    
    CALL deduce_optimal_histogram_properties
    
    OPEN (10, FILE=filename, FORM='formatted')
    WRITE (10,*) '# r  rdf'

    error%found = .FALSE.    
    IF (nb_species/=1) STOP "compute_rdf.f90 is written for 1 solvent species only."

    ALLOCATE (rdf (nbins) ,SOURCE=0._dp)
    
    ! CALL UTest_histogram_3D
    DO n = 1, SIZE(soluteSite) ! loop over all sites of the solute
        CALL histogram_3d (array(:,:,:), rdf)
        IF (error%found) THEN; PRINT*,error%msg; STOP; END IF
        ! Write to output/rdf.out
        WRITE(10,*)'# solute site', n
        WRITE(10,*) 0._dp,0._dp ! we impose
        DO bin =1,nbins
            WRITE(10,*) (REAL(bin,dp)-0.5_dp)*dr, rdf(bin) ! For bin that covers 0<r<dr, I print at 0.5dr, i.e., at the middle of the bin
        END DO
        WRITE(10,*)
    END DO
    DEALLOCATE (rdf)
    CLOSE(10)

    
    CONTAINS

        !===========================================================================================================================
        SUBROUTINE deduce_optimal_histogram_properties
        !===========================================================================================================================
            IMPLICIT NONE
            nbins = CEILING (2._sp*REAL (PRODUCT (spaceGrid%n_nodes), sp)**(1._sp/3._sp)) ! Rice Rule, see http://en.wikipedia.org/wiki/Histogram
            RdfMaxRange= SQRT(3._dp)*MINVAL(spaceGrid%length)/2._dp
            dr = RdfMaxRange/REAL(nbins,dp) ! Width of each bin of the histogram
        END SUBROUTINE deduce_optimal_histogram_properties
        !===========================================================================================================================
    
        !===========================================================================================================================
        SUBROUTINE histogram_3d (array3D, rdf)
        !===========================================================================================================================
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: array3D(spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3))
            REAL(dp), INTENT(OUT) :: rdf(nbins)
            INTEGER(i2b), ALLOCATABLE :: recurrence_bin(:)
            REAL(dp) :: r
            INTEGER(i2b) :: bin, i, j, k
            
            ALLOCATE( recurrence_bin (nbins) ,SOURCE=0)
            rdf = 0._dp
            
            ! Transform array(position) in rdf(radialdistance)
            ! counts the total number of appearence of a value in each bin
            DO CONCURRENT (i=1:spaceGrid%n_nodes(1), j=1:spaceGrid%n_nodes(2), k=1:spaceGrid%n_nodes(3))
                r = NORM2(REAL([i,j,k]-1,dp)*spaceGrid%dl - soluteSite(n)%r)
                bin = INT(r/dr)+1
                IF (bin>nbins) THEN
                    CYCLE
                ELSE IF (bin<=0) THEN
                    error%found = .TRUE.
                    error%msg = "I found a null or negative bin index in compute_rdf.f90.\\ STOP"
                END IF
                recurrence_bin(bin) =recurrence_bin(bin) +1
                rdf           (bin) =rdf(bin) +array3D(i,j,k)
            END DO

            ! test if a problem found
            IF (error%found) THEN
                PRINT*, error%msg
                STOP
            END IF

            ! normalize the rdf
            WHERE (recurrence_bin/=0)
                rdf = rdf/REAL(recurrence_bin,dp)
            ELSEWHERE
                rdf = 0.0_dp
            END WHERE
            
            DEALLOCATE (recurrence_bin)
        END SUBROUTINE
        !===========================================================================================================================

        !===========================================================================================================================
        SUBROUTINE UTest_histogram_3D
            IMPLICIT NONE
            REAL(dp), ALLOCATABLE :: nullarray3D (:,:,:) 
            REAL(dp), ALLOCATABLE :: rdfUT(:)
            INTEGER(i2b) :: bin

            ALLOCATE (nullarray3D (spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3)) ,SOURCE=0._dp)
            ALLOCATE (rdfUT(nbins),source=0._dp)

            ! Test 1
            nullarray3D = 0._dp
            CALL histogram_3d (nullarray3D, rdfUT)
            IF (ANY(rdfUT/=0._dp)) STOP "Test 1 in UTest_histogram_3D not passed."

            ! Test 2
            nullarray3D = 100._dp
            CALL histogram_3d (nullarray3D, rdfUT)
            IF (ANY(rdfUT<0._dp)) STOP "Test 2 in UTest_histogram_3D not passed."

            ! Test 3
            nullarray3D = -1._dp
            CALL histogram_3d (nullarray3D, rdfUT)
            IF (ANY(rdfUT>0._dp)) STOP "Test 3 in UTest_histogram_3D not passed."

            ! Test 4
            nullarray3D = 100._dp
            CALL histogram_3d (nullarray3D, rdfUT)
            IF (ANY(rdfUT>100._dp)) STOP "Test 4 in UTest_histogram_3D not passed."

            DEALLOCATE (nullarray3D, rdfUT)
        END SUBROUTINE UTest_histogram_3D
        !===========================================================================================================================

END SUBROUTINE compute_rdf
