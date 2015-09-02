! Given a 3 dimensional scalar field, this subroutine prints the value of this field in a given plane.
! For now, only planes perpendicular to the x, y and z directions are allowed.
! TODO: should be extended to any orientation of the plane.

SUBROUTINE compute_planar_density(array,filename)

    use precision_kinds ,ONLY: dp,i2b
    use system          ,ONLY: solute
    use module_grid, only: grid

    IMPLICIT NONE

    CHARACTER(50), INTENT(IN)   :: filename
    REAL(dp), INTENT(IN)        :: array(grid%n_nodes(1),grid%n_nodes(2),grid%n_nodes(3))
    INTEGER(i2b)                :: plandir,ig,i,j,k
    REAL(dp)                    :: x_com(grid%n_nodes(1)), y_com(grid%n_nodes(2)), z_com(grid%n_nodes(3))

    INTEGER(i2b) :: nfft1, nfft2, nfft3
    REAL(dp) :: lx, ly, lz

    lx= grid%length(1)
    ly= grid%length(2)
    lz= grid%length(3)
    nfft1= grid%n_nodes(1)
    nfft2= grid%n_nodes(2)
    nfft3= grid%n_nodes(3)


    ! For now that only planes perpendicular to x, y and z are  Identify the plan coordinate which is 0 (see restrictions to using this program for now)
    IF      (solute%site(1)%r(1)==solute%site(2)%r(1) .and. solute%site(1)%r(1)==solute%site(3)%r(1)) THEN
        plandir=1
    ELSE IF (solute%site(1)%r(2)==solute%site(2)%r(2) .and. solute%site(1)%r(2)==solute%site(3)%r(2)) THEN
        plandir=2
    ELSE IF (solute%site(1)%r(3)==solute%site(2)%r(3) .and. solute%site(1)%r(3)==solute%site(3)%r(3)) THEN
        plandir=3
    ELSE
        RETURN
    END IF

    ! Get its grid index called id
    IF (plandir==1) THEN
        ig= NINT(solute%site(1)%r(1)*REAL(nfft1,dp)/Lx) +1
    ELSE IF (plandir==2) then
        ig= NINT(solute%site(1)%r(2)*REAL(nfft2,dp)/Ly) +1
    ELSE IF (plandir==3) then
        ig= NINT(solute%site(1)%r(3)*REAL(nfft3,dp)/Lz) +1
    ELSE
        print*,'I did not find the direction of the plane in compute_planar_density.f90'
        return
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
                    WRITE(10,100) y_com(j), z_com(k), array(ig,j,k)
                END DO
            END DO
        CASE (2)
            DO i=1,nfft1
                DO k=1,nfft3
                    WRITE(10,100)x_com(i),z_com(k),array(i,ig,k)
                END DO
            END DO
        CASE (3)
            DO i=1,nfft1
                DO j=1,nfft2
                    WRITE(10,100)x_com(i),y_com(j),array(i,j,ig)
                END DO
            END DO
        END SELECT
    CLOSE(10)

END SUBROUTINE
