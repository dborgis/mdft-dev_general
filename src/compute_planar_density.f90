! Given a 3 dimensional scalar field, this subroutine prints the value of this field in a given plane.
! For now, only planes perpendicular to the x, y and z directions are allowed.
! TODO: should be extended to any orientation of the plane.

SUBROUTINE compute_planar_density(array,filename)

    USE precision_kinds ,ONLY: dp,i2b
    USE system          ,ONLY: nfft1,nfft2,nfft3,nb_solute_sites,Lx,Ly,Lz,spaceGrid, soluteSite

    IMPLICIT NONE
    
    CHARACTER(50), INTENT(IN)   :: filename
    REAL(dp), INTENT(IN)        :: array(nfft1,nfft2,nfft3)
    INTEGER(i2b)                :: plandir,id,i,j,k
    REAL(dp)                    :: x_com(nfft1), y_com(nfft2), z_com(nfft3)
    
    ! For now that only planes perpendicular to x, y and z are  Identify the plan coordinate which is 0 (see restrictions to using this program for now)
    IF      (soluteSite(1)%r(1)==soluteSite(2)%r(1) .and. soluteSite(1)%r(1)==soluteSite(3)%r(1)) THEN
        plandir=1
    ELSE IF (soluteSite(1)%r(2)==soluteSite(2)%r(2) .and. soluteSite(1)%r(2)==soluteSite(3)%r(2)) THEN
        plandir=2
    ELSE IF (soluteSite(1)%r(3)==soluteSite(2)%r(3) .and. soluteSite(1)%r(3)==soluteSite(3)%r(3)) THEN
        plandir=3
    ELSE
        RETURN
    END IF

    ! Get its grid index called id
    IF (plandir==1) THEN
        id= NINT(soluteSite(1)%r(1)*REAL(nfft1,dp)/Lx) +1
    ELSE IF (plandir==2) then
        id= NINT(soluteSite(1)%r(2)*REAL(nfft2,dp)/Ly) +1
    ELSE IF (plandir==3) then
        id= NINT(soluteSite(1)%r(3)*REAL(nfft3,dp)/Lz) +1
    ELSE
        WRITE(*,*) 'I did not find the direction of the plane in compute_planar_density.f90'
    END IF

    ! Compute grid points cartesian coordinates
    FORALL(i=1:nfft1) x_com(i)=REAL(i-1,dp)*Lx/REAL(nfft1,dp)
    FORALL(j=1:nfft2) y_com(j)=REAL(j-1,dp)*Ly/REAL(nfft2,dp)
    FORALL(k=1:nfft3) z_com(k)=REAL(k-1,dp)*Lz/REAL(nfft3,dp)
    
    ! Print density in this plan
    OPEN(10,FILE=filename,FORM='formatted')
        100 FORMAT (3(F10.5))
        SELECT CASE (plandir)
        CASE (1)
            DO j=1,nfft2
                DO k=1,nfft3
                    WRITE(10,100) y_com(j), z_com(k), array(id,j,k)
                END DO
            END DO
        CASE (2)
            DO i=1,nfft1
                DO k=1,nfft3
                    WRITE(10,100)x_com(i),z_com(k),array(i,id,k)
                END DO
            END DO
        CASE (3)
            DO i=1,nfft1
                DO j=1,nfft2
                    WRITE(10,100)x_com(i),y_com(j),array(i,j,id)
                END DO
            END DO
        END SELECT
    CLOSE(10)

END SUBROUTINE
