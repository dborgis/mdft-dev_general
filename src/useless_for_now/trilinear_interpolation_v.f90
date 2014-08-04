! on a un champ V(nfft1, nfft2, nfft3) sur une grille de taille nfft1 x nfft2 x nfft3
! Connaissant ce champ en des points discrets, on veut le calculer en un point quelconque de coordonnees x, y, z
! C'est ce que fait cette routine par interpolation trilineaire
! La stratégie est la suivante :
! 1/ on traduit les coordonnées cartesiennes x y z en coordonnées dans l'espace des indices (xi, yi, zi)
! 2/ on détermine les huit points de la grille discrete qui entourent notre point
! 3/ on en déduit l'interpolation trilineaire. La forme très simple utilisée pour le calcul est due au travail dans l'espace des indices
! ainsi on a toujours un pas entre les indices discrets = 1.
SUBROUTINE trilinear_interpolation_v(x, y, z, V_interpolated)
 USE precision_kinds,only : dp , i2b
 use system, only : V_coulomb, nfft1, nfft2, nfft3, Lx, Ly, Lz, deltax, deltay, deltaz
 IMPLICIT NONE
! real(dp), dimension (nfft1, nfft2, nfft3), intent(in) :: V_grid ! potentiel connu sur une grille discrete
 real(dp), intent(in) :: x,y,z ! coordonnee cartesienne du point dont on cherche V par interpolation trilineaire
 real(dp), intent(out) :: V_interpolated ! potentiel interpolé au point qui nous interesse de coordonnees cartesiennes x y z
 real(dp) :: xi, yi, zi ! indice (réel puisqu'on est en fait entre deux indices) correspondant a x
 ! floor(xic) est le kind real de int(xic) == real(int(xi))point de la grille juste en dessous de x : rappel que int(9.9) = 9
 ! point de la grille juste en dessous de y
 ! point de la grille juste en dessous de z
 real(dp) :: t, u, v ! qui devrait etre x mod(1) == mod(x,1) en fortran a verifier plus tard
 real(dp) :: omt, omu, omv ! 1-t
 ! Vi potentiel en chacun des huits points qui entourent notre point de coordonnees x y z
 real(dp) :: V1, V2, V3, V4, V5, V6, V7, V8
 integer(i2b) :: j, k, l ! dummy index
 integer(i2b) :: jp1,kp1,lp1 ! j+1
 
 xi = modulo(1.0_dp + x/deltax, Lx) ! indice (réel puisqu'on est en fait entre 8 pts de grille) correspondant a x
 yi = modulo(1.0_dp + y/deltay, Ly)
 zi = modulo(1.0_dp + z/deltaz, Lz)
 ! floor(xic) est le kind real de int(xic) == real(int(xi))point de la grille juste en dessous de x : rappel que int(9.9) = 9
 ! point de la grille juste en dessous de y
 ! point de la grille juste en dessous de z
 t = (xi-aint(xi)) ! qui devrait etre x mod(1) == mod(x,1) en fortran a verifier plus tard
 u = (yi-aint(yi)) ! ANCIENNEMENT FLOOR
 v = (zi-aint(zi))
 omt = 1.0_dp-t
 omu = 1.0_dp-u
 omv = 1.0_dp-v
! 8 points qui entourent le point dont on veut l'interpolation
 j = int(xi)
 k = int(yi)
 l = int(zi)
 if (j==0) then
  j=nfft1
  jp1=1
 ELSE IF (j>0 .and. j<nfft1) then
  jp1 = j+1
 ELSE IF (j==nfft1) then
  jp1 = 1
 ELSE IF (j>nfft1) then
  j=j-nfft1
  jp1 = j+1
 END IF
 if (k==0) then
  k=nfft2
  kp1=1
 ELSE IF (k>0 .and. k<nfft2) then
  kp1 = k+1
 ELSE IF (k==nfft2) then
  kp1 = 1
 ELSE IF (k>nfft2) then
  k=k-nfft2
  kp1 = k+1
 END IF
 if (l==0) then
  l=nfft3
  lp1=1
 ELSE IF (l>0 .and. l<nfft3) then
  lp1 = l+1
 ELSE IF (l==nfft3) then
  lp1 = 1
 ELSE IF (l>nfft3) then
  l=l-nfft3
  lp1 = l+1
 END IF
 ! Vi potentiel en chacun des huits points qui entourent notre point de coordonnees x y z
 V1 = V_coulomb(j,k,l)
 V2 = V_coulomb(jp1,k,l)
 V3 = V_coulomb(j,kp1,l)
 V4 = V_coulomb(j,k,lp1)
 V5 = V_coulomb(jp1,kp1,l)
 V6 = V_coulomb(j,kp1,lp1)
 V7 = V_coulomb(jp1,k,lp1)
 V8 = V_coulomb(jp1,kp1,lp1)
 
 V_interpolated = omt*omu*omv*V1 + t*omu*omv*V2 + omt*u*omv*V3 + omt*omu*v*V4 + t*u*omv*V5 + omt*u*v*V6 + t*omu*v*V7 + t*u*v*V8
END SUBROUTINE trilinear_interpolation_v
