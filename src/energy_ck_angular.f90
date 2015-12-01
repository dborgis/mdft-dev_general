SUBROUTINE energy_ck_angular (Fexc_ck_angular)
!
    use precision_kinds, ONLY : i2b, dp
!     use system,          ONLY : thermocond, grid
!       use module_solvent, only: solvent
!     use quadrature,      ONLY : Omx, Omy, Omz, angGrid, molRotGrid, molRotSymOrder
!     use module_minimizer,       ONLY : cg_vect_new, FF, dF_new
!     use constants,       ONLY : twopi
!     use fft,             ONLY : fftw3, kx, ky, kz
!     use module_input,           ONLY : getinput%log
!     use dcf,             ONLY : ck_angular, angleInd, angleVal, delta_k, delta_k_ck_angular, num_phi, num_cos, num_psi, &
!                                 c_s, c_d, c_delta, c_q, nb_k
!
    IMPLICIT NONE
!
    REAL(dp),  INTENT(OUT) :: Fexc_ck_angular ! excess part of the free energy
!     REAL(dp),    PARAMETER :: root3=SQRT(3._dp),root03=SQRT(0.3_dp)
!     COMPLEX(dp), PARAMETER :: imag=(0._dp,1._dp)
!     INTEGER(i2b) :: icg, i, j, k, l, m, n, o, p, o1, p1, o2, p2, nfft1, nfft2, nfft3, s
!     INTEGER(i2b) :: ik, iphi, ipsi1, ipsi2, icos1, icos2
!     INTEGER(i2b) :: n1, n2, n3, n4, n5
!     INTEGER(i2b) :: karim
!     REAL(dp) :: chi
!     REAL(dp) :: Nk, integrationFactor, local_weight
!     REAL(dp) :: om1_dot_om2, k_dot_om1, k_dot_om2
!     REAL(dp) :: phi12_value, psi1_value, psi2_value, cos1_value, cos2_value
!     REAL(dp) :: delta_phi, delta_psi, delta_cos
!     REAL(dp), DIMENSION(0:1) :: deltaphi, deltapsi1, deltapsi2, deltacos1, deltacos2
!     REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: delta_rho, gamma
!     COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: delta_rho_k
!     COMPLEX(dp) :: Ck
!
! ! Initialisation
!     nfft1 = grid%n_nodes(1); nfft2 = grid%n_nodes(2); nfft3 = grid%n_nodes(3)
!     Nk = REAL(nfft1*nfft2*nfft3,dp) ! Total number of k grid points in real
!     Fexc_ck_angular = 0.0_dp
!     integrationFactor = grid%dv * solvent(1)%rho0
!     delta_phi = twopi / num_phi
!     delta_psi = twopi / (num_psi * molRotSymOrder)
!     delta_cos = 2._dp / num_cos
!
!     ALLOCATE(delta_rho  (nfft1    , nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles), SOURCE=0._dp        )
!     ALLOCATE(delta_rho_k(nfft1/2+1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles), SOURCE=(0._dp,0._dp))
!
!     delta_rho(i,j,k,o,p) = cg_vect_new(i,j,k,o,p,1)**2 -1._dp
!
! ! Fourier transform of delta_rho to delta_rho_k
!     DO o1 = 1, angGrid%n_angles; DO p1 = 1, molRotGrid%n_angles
!         fftw3%in_forward = delta_rho(:,:,:,o1,p1)
!         CALL dfftw_execute(fftw3%plan_forward)
!         delta_rho_k(:,:,:,o1,p1) = fftw3%out_forward
!     END DO; END DO
!
!     DEALLOCATE(delta_rho)
!     ALLOCATE(gamma(nfft1, nfft2, nfft3, angGrid%n_angles, molRotGrid%n_angles), SOURCE=0._dp)
!
! ! Logical for choice of method
!     karim = 0
!     IF (getinput%log("ck_angular"))        karim = 1
!     IF (getinput%log("ck_debug")) THEN;    karim = 2; delta_k_ck_angular = delta_k; END IF
!     IF (getinput%log("ck_debug_extended")) karim = 3
!     IF (getinput%log("ck_angular") .AND. getinput%log("ck_angular_interpolation")) karim = 4
!     IF (karim == 0) STOP 'FATAL ERROR: Karim est nul in energy_ck_angular 0_0'
!
! ! Cauculation of gamma(r,Omega)
!     DO o1 = 1, angGrid%n_angles; DO p1 = 1, molRotGrid%n_angles ! not DO CONCURRENT because of fftw3
!         fftw3%in_backward = (0._dp,0._dp) ! fftw3%in_backward is gamma_k(:,:,:,o1,p1)
!
!         DO CONCURRENT (l=1:nfft1/2+1, m=1:nfft2, n=1:nfft3, o2=1:angGrid%n_angles, p2=1:molRotGrid%n_angles)
!             ik = INT((kx(l)**2 + ky(m)**2 + kz(n)**2)**0.5_dp / delta_k_ck_angular + 0.5_dp) + 1 ! index k
!             IF(ik > nb_k) karim = 0
!
!             SELECT CASE(karim)
!             CASE(0) ! ik exceeds the maximum
!                 PRINT*, 'WARNING: Index of k ',ik,' exceeds its maximum ',nb_k,', ck_angular is put to zero.'
!                 ck = (0._dp,0._dp)
!
!             CASE(1,4) ! ck_angular, ck_angular_interpolation
!                 IF (karim==4) THEN ! linear interpolation
!                     phi12_value = MODULO(angleVal(l,m,n,o1,p1)%phi - angleVal(l,m,n,o2,p2)%phi,twopi) ! = phi1 - phi2
!                     psi1_value  = angleVal(l,m,n,o1,p1)%psi
!                     psi2_value  = angleVal(l,m,n,o2,p2)%psi
!                     cos1_value  = angleVal(l,m,n,o1,p1)%costheta
!                     cos2_value  = angleVal(l,m,n,o2,p2)%costheta
!
!                     iphi = INT(phi12_value/delta_phi + 0.5_dp) ! Interpolation always pick the clockwise nearest angle
!                     ipsi1 = INT(psi1_value/delta_psi + 0.5_dp)
!                     ipsi2 = INT(psi2_value/delta_psi + 0.5_dp)
!                     icos1 = INT((cos1_value + 1._dp)/delta_cos + 0.5_dp) ! When icos == 0 or num_cos - 1, extrapolation needed
!                     icos2 = INT((cos2_value + 1._dp)/delta_cos + 0.5_dp)
!
!                     deltaphi(1)  = (phi12_value - (iphi - 0.5_dp) * delta_phi) / delta_phi ! Calculation of weights
!                     deltaphi(0)  = 1._dp - deltaphi(1)
!                     deltapsi1(1) = (psi1_value - (ipsi1 - 0.5_dp) * delta_psi) / delta_psi
!                     deltapsi1(0) = 1._dp - deltapsi1(1)
!                     deltapsi2(1) = (psi2_value - (ipsi2 - 0.5_dp) * delta_psi) / delta_psi
!                     deltapsi2(0) = 1._dp - deltapsi2(1)
!
!                     iphi  = MODULO(iphi  - 1, num_phi) + 1 ! To avoid angle index 0
!                     ipsi1 = MODULO(ipsi1 - 1, num_psi) + 1
!                     ipsi2 = MODULO(ipsi2 - 1, num_psi) + 1
!
!                     IF (icos1 == 0) icos1 = 1 ! Linear extrapolation for cos_values: same formula and sample points as interpolation.
!                     IF (icos2 == 0) icos2 = 1
!                     IF (icos1 == num_cos) icos1 = num_cos - 1
!                     IF (icos2 == num_cos) icos2 = num_cos - 1
!
!                     deltacos1(1) = (cos1_value + 1._dp - (icos1 - 0.5_dp) * delta_cos) / delta_cos
!                     deltacos1(0) = 1._dp - deltacos1(1)
!                     deltacos2(1) = (cos2_value + 1._dp - (icos2 - 0.5_dp) * delta_cos) / delta_cos
!                     deltacos2(0) = 1._dp - deltacos2(1)
!
!                     Ck = 0._dp
!                     DO n1=0,1; DO n2=0,1; DO n3=0,1; DO n4=0,1; DO n5=0,1
!                         Ck = Ck + deltapsi1(n1) * deltapsi2(n2) * deltaphi(n3) * deltacos1(n4) * deltacos2(n5) &
!                            * ck_angular(MOD(ipsi1+n1-1,num_psi)+1, MOD(ipsi2+n2-1,num_psi)+1, MOD(iphi+n3-1,num_phi)+1, &
!                              icos1+n4, icos2+n5, ik)
!                     END DO; END DO; END DO; END DO; END DO
!
!                 ELSE ! zero order interpolation
!                     phi12_value = MODULO(angleInd(l,m,n,o1,p1)%phi - angleInd(l,m,n,o2,p2)%phi,twopi) ! Luc Phi1 - Phi2, Daniel Phi2 - Phi1
!                     iphi = MOD(INT(phi12_value*num_phi/twopi), num_phi) + 1
!                     ipsi1 = angleInd(l,m,n,o1,p1)%psi
!                     ipsi2 = angleInd(l,m,n,o2,p2)%psi
!                     icos1 = angleInd(l,m,n,o1,p1)%costheta
!                     icos2 = angleInd(l,m,n,o2,p2)%costheta
!
!                     Ck = ck_angular(ipsi1,ipsi2,iphi,icos1,icos2,ik)
!                 END IF
!
!             ! CASE(2) ! ck_debug
!             !   stop "maybe ck_debug doesnt work anymore since c_s c_delta and c_d have been improved a lot in 2015"
!             !     om1_dot_om2 = Omx(o1)*Omx(o2) + Omy(o1)*Omy(o2) + Omz(o1)*Omz(o2)
!             !     IF ((kx(l)**2 + ky(m)**2 + kz(n)**2)**0.5_dp == 0._dp) THEN
!             !         k_dot_om1 = Omz(o1)
!             !         k_dot_om2 = Omz(o2)
!             !     ELSE
!             !         k_dot_om1 = (kx(l)*Omx(o1) + ky(m)*Omy(o1) + kz(n)*Omz(o1))/(kx(l)**2 + ky(m)**2 + kz(n)**2)**0.5_dp
!             !         k_dot_om2 = (kx(l)*Omx(o2) + ky(m)*Omy(o2) + kz(n)*Omz(o2))/(kx(l)**2 + ky(m)**2 + kz(n)**2)**0.5_dp
!             !     END IF
!             !
!             !     Ck = c_s(ik) + c_delta(ik) * om1_dot_om2 + c_d(ik) * (3._dp * k_dot_om1 * k_dot_om2 - om1_dot_om2) ! For factors of Daniel
!
!             ! CASE(3) ! ck_debug_extended
!             !   stop "maybe ck_debug does not work anymore since c_s c_delta and c_d have been improved a lot since dec 2014"
!             !
!             !     om1_dot_om2 = Omx(o1)*Omx(o2) + Omy(o1)*Omy(o2) + Omz(o1)*Omz(o2)
!             !     IF ((kx(l)**2 + ky(m)**2 + kz(n)**2)**0.5_dp == 0._dp) THEN
!             !         k_dot_om1 = Omz(o1)
!             !         k_dot_om2 = Omz(o2)
!             !     ELSE
!             !         k_dot_om1 = (kx(l)*Omx(o1) + ky(m)*Omy(o1) + kz(n)*Omz(o1))/(kx(l)**2 + ky(m)**2 + kz(n)**2)**0.5_dp
!             !         k_dot_om2 = (kx(l)*Omx(o2) + ky(m)*Omy(o2) + kz(n)*Omz(o2))/(kx(l)**2 + ky(m)**2 + kz(n)**2)**0.5_dp
!             !     END IF
!             !
!             !     Ck = c_s(ik) - root3 * c_delta(ik) * om1_dot_om2 &
!             !         + imag * c_q(ik) * k_dot_om1 - imag * c_q(ik) * k_dot_om2 & ! here ck is not the conjugated function of Luc
!             !         + c_d(ik) * root03 * (3._dp * k_dot_om1 * k_dot_om2 - om1_dot_om2) ! For factors of Luc
!
!             END SELECT
!
!             fftw3%in_backward(l,m,n) = fftw3%in_backward(l,m,n) &
!                                      + Ck * delta_rho_k(l,m,n,o2,p2) * angGrid%weight(o2) * molRotGrid%weight(p2)
!         END DO
!
!     ! Transform gamma(k,Omega) to gamma(r,Omega)
!         CALL dfftw_execute(fftw3%plan_backward)
!         gamma(:,:,:,o1,p1) = fftw3%out_backward/Nk *(-thermocond%kbT)* solvent(1)%rho0
!
!     END DO; END DO
!
! ! Calculation of excess part of the free energy functional
!     icg = 0
!     s = 1
!     DO i = 1, nfft1; DO j = 1, nfft2; DO k = 1, nfft3; DO o = 1, angGrid%n_angles; DO p = 1, molRotGrid%n_angles
!         icg = icg + 1
!         chi = cg_vect_new(i,j,k,o,p,1)
!         local_weight = angGrid%weight(o) * molRotGrid%weight(p)
!         Fexc_ck_angular = Fexc_ck_angular + gamma(i,j,k,o,p) * (chi ** 2 - 1) * local_weight
!         dF_new(i,j,k,o,p,s) = dF_new(i,j,k,o,p,s) + 2._dp * chi * gamma(i,j,k,o,p) * integrationFactor * local_weight
!     END DO; END DO; END DO; END DO; END DO
!
!     Fexc_ck_angular = Fexc_ck_angular / 2._dp * integrationFactor
!     FF = FF + Fexc_ck_angular
!
!     DEALLOCATE(delta_rho_k); DEALLOCATE(gamma)
!
! END SUBROUTINE energy_ck_angular
!
! !-----------------------------------------------------------------------------------------------------------------------------------
!
! SUBROUTINE energy_ck_en_harm_sph! (Fexc_ck_angular)
!
!     STOP 'For computing energy by expansion on spherical harmonics with Luc''s full projections, to be developed.'
!
! END SUBROUTINE energy_ck_en_harm_sph
!
! !-----------------------------------------------------------------------------------------------------------------------------------
!
! SUBROUTINE proj_en_harm_sph!(delta_rho_k, delta_rho_proj)
!
!     STOP 'For calculating projections delta_rho_proj(k,m,mu,mu'') from delta_rho_k(k,omega), to be developed.'
!
! END SUBROUTINE proj_en_harm_sph
!
! !-----------------------------------------------------------------------------------------------------------------------------------
!
! SUBROUTINE rot_proj!(proj, proj_loc)
!
!     STOP 'For rotating projections to local coordinates or inverse, to be developed.'
!
! END SUBROUTINE rot_proj
END SUBROUTINE energy_ck_angular
