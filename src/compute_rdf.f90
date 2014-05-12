! I compute the radial distribution function (rdf) around each solute site
SUBROUTINE compute_rdf (array,filename)

    USE precision_kinds ,ONLY: dp,i2b
    USE system          ,ONLY: x_mol,y_mol,z_mol,id_mol,id_solv,nb_species,spaceGrid,soluteSite
    USE input           ,ONLY: input_dp,input_int
    
    IMPLICIT NONE

    REAL(dp), DIMENSION(spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3),nb_species), INTENT(IN) :: array
    CHARACTER( 50 ) , intent(in) :: filename
    INTEGER(i2b):: nb_id_solv,nb_id_mol,nbins
    REAL(dp):: RdfMaxRange,dr,r,xnm2,ynm2,znm2,dx,dy,dz
    REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: rdf
    INTEGER(i2b), ALLOCATABLE, DIMENSION (:,:) :: recurrence_bin
    INTEGER(i2b):: i,j,k,n,bin,s,nfft1,nfft2,nfft3
    LOGICAL :: foundErr

    nb_id_solv = MAXVAL(id_solv)
    nb_id_mol  = MAXVAL(id_mol)
    nbins=input_int('nbinsrdf') ! The rdf is a histogram of this number of bins
    RdfMaxRange=input_dp('rdfmaxrange') ! max range of the rdf (10 Ang is a minimum)
    !vRdfMaxRange=min(Lx,Ly,Lz)*sqrt(3.0_dp) !> Maximum distance to which calculate radial distribution function
    dr = RdfMaxRange/REAL(nbins,dp) ! Width of each bin of the histogram

    OPEN(10,FILE=filename,FORM='formatted')
    
    100 FORMAT(2(F10.5))
    WRITE(10,*)'# r  rdf'

    ALLOCATE( recurrence_bin ( SIZE(soluteSite), 1:nbins) ,SOURCE=0)
    ALLOCATE( rdf            ( SIZE(soluteSite), 1:nbins) ,SOURCE=0._dp)

    nfft1 =spaceGrid%n_nodes(1)
    nfft2 =spaceGrid%n_nodes(2)
    nfft3 =spaceGrid%n_nodes(3)
    dx =spaceGrid%dl(1)
    dy =spaceGrid%dl(2)
    dz =spaceGrid%dl(3)
    foundErr=.FALSE.
    
    DO s=1,nb_species
        IF (nb_species/=1) STOP "compute_rdf.f90 valid only for 1 solvent molecule"
        recurrence_bin = 0
        rdf = 0.0_dp
        WRITE(10,*)'# You want',s,'molecular solvents'
        
        ! Transform array(position) in rdf(radialdistance)
        ! counts the total number of appearence of a value in each bin
        DO CONCURRENT (n=1:SIZE(soluteSite), i=1:nfft1, j=1:nfft2, k=1:nfft3)
            xnm2 =((i-1)*dx-x_mol(n))**2
            ynm2 =((j-1)*dy-y_mol(n))**2
            znm2 =((k-1)*dz-z_mol(n))**2
            r    =SQRT(xnm2+ynm2+znm2)
            bin  =INT(r/dr)+1
            IF (bin>nbins) THEN
                CYCLE
            ELSE IF (bin<=0) THEN
                foundErr =.TRUE.
            END IF
            recurrence_bin(n,bin) =recurrence_bin(n,bin) +1
            rdf           (n,bin) =rdf(n,bin) +array(i,j,k,s)
        END DO
        
        IF (foundErr) STOP "I found a null or negative bin index in compute_rdf.f90.\\ STOP"

        !normalize the rdf
        WHERE (recurrence_bin/=0)
            rdf = rdf/REAL(recurrence_bin,dp)
        ELSEWHERE
            rdf = 0.0_dp
        END WHERE
    
        ! Write rdf
        DO n=1,SIZE(soluteSite)
            WRITE(10,*)'#site ' , n
            WRITE(10,100)0._dp,0._dp
            DO bin =1,nbins
                WRITE(10,100) (REAL(bin,dp)-0.5_dp)*dr, rdf(n,bin) ! For bin that covers 0<r<dr, I print at 0.5dr, i.e., at the middle of the bin
            END DO
            WRITE(10,*)
        END DO
        WRITE(10,*)! white line between solvent molecules
    END DO
    
    CLOSE(10)

END SUBROUTINE compute_rdf
