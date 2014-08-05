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

    REAL(dp), DIMENSION(spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3),nb_species), INTENT(IN) :: array
    CHARACTER(50), INTENT(IN) :: filename
    REAL(dp) :: RdfMaxRange, dr, r
    REAL(dp), ALLOCATABLE :: rdf(:,:)
    INTEGER(i2b), ALLOCATABLE :: recurrence_bin(:,:)
    INTEGER(i2b):: i,j,k,n,bin,s,nbins
    TYPE :: errortype
        LOGICAL :: found
        CHARACTER(180) :: msg
    END TYPE
    TYPE (errortype) :: error
    
    nbins = CEILING (2._sp*REAL (PRODUCT (spaceGrid%n_nodes), sp)**(1._sp/3._sp)) ! Rice Rule, see http://en.wikipedia.org/wiki/Histogram
    RdfMaxRange= SQRT(3._dp)*MINVAL(spaceGrid%length)/2._dp
    dr = RdfMaxRange/REAL(nbins,dp) ! Width of each bin of the histogram

    OPEN(10,FILE=filename,FORM='formatted')
        WRITE(10,*)'# r  rdf'

    ALLOCATE( recurrence_bin ( SIZE(soluteSite), nbins) ,SOURCE=0)
    ALLOCATE( rdf            ( SIZE(soluteSite), nbins) ,SOURCE=0._dp)

    error%found = .FALSE.
    
    DO s=1,nb_species
        IF (nb_species/=1) STOP "compute_rdf.f90 is written for 1 solvent species only."
        recurrence_bin = 0
        rdf = 0.0_dp
        WRITE(10,*)'# You want',s,'molecular solvents'

        ! Transform array(position) in rdf(radialdistance)
        ! counts the total number of appearence of a value in each bin
        DO CONCURRENT (n=1:SIZE(soluteSite), i=1:spaceGrid%n_nodes(1), j=1:spaceGrid%n_nodes(2), k=1:spaceGrid%n_nodes(3))
            r = NORM2(REAL([i,j,k]-1,dp)*spaceGrid%dl - soluteSite(n)%r)
            bin = INT(r/dr)+1
            IF (bin>nbins) THEN
                CYCLE
            ELSE IF (bin<=0) THEN
                error%found = .TRUE.
                error%msg = "I found a null or negative bin index in compute_rdf.f90.\\ STOP"
            END IF
            recurrence_bin(n,bin) =recurrence_bin(n,bin) +1
            rdf           (n,bin) =rdf(n,bin) +array(i,j,k,s)
        END DO

        IF (error%found) THEN
            PRINT*,error%msg
            STOP
        END IF

        !normalize the rdf
        WHERE (recurrence_bin/=0)
            rdf = rdf/REAL(recurrence_bin,dp)
        ELSEWHERE
            rdf = 0.0_dp
        END WHERE
    
        ! Write rdf
        DO n=1,SIZE(soluteSite)
            WRITE(10,*)'#site ' , n
            WRITE(10,*)0._dp,0._dp
            DO bin =1,nbins
                WRITE(10,*) (REAL(bin,dp)-0.5_dp)*dr, rdf(n,bin) ! For bin that covers 0<r<dr, I print at 0.5dr, i.e., at the middle of the bin
            END DO
            WRITE(10,*)
        END DO
        WRITE(10,*)! white line between solvent molecules
    END DO
    
    DEALLOCATE (recurrence_bin, rdf)
    CLOSE(10)

END SUBROUTINE compute_rdf
