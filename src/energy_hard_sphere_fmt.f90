SUBROUTINE energy_hard_sphere_fmt (Ffmt)

    USE precision_kinds  ,ONLY: dp , i2b
    USE system           ,ONLY: thermocond, nb_species, mole_fraction, spaceGrid, solvent
    USE quadrature       ,ONLY: molRotSymOrder, angGrid, molRotGrid
    USE minimizer        ,ONLY: cg_vect , FF , dF
    USE constants        ,ONLY: pi , FourPi , twopi, zeroC
    USE fft              ,ONLY: fftw3
    USE input            ,ONLY: input_char
    USE hardspheres      ,ONLY: hs

    IMPLICIT NONE

    REAL(dp), INTENT(OUT) :: Ffmt ! Internal part of free energy
    INTEGER(i2b) :: icg,i,j,k,o,p,s,nfft1,nfft2,nfft3,wdl,wdu
    REAL(dp) :: local_density, psi, dV
    REAL(dp) :: nb_molecules(size(solvent))
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: rho ! density per angle (recall : rho_0 = n_0 / 4pi ) ! x y z nb_species
    REAL(dp) :: time0 , time1 ! time stamps
    REAL(dp) :: w0,w1,w2,w3,omw3, kT
    REAL(dp)   , ALLOCATABLE , DIMENSION(:,:,:,:) :: wd ! weighted density at node i,j,k, for index 0:3
    REAL(dp)   , ALLOCATABLE , DIMENSION(:,:,:) :: one_min_wd_3 ! dummy for 1 - weighted density 3
    REAL(dp)   , ALLOCATABLE , DIMENSION(:,:,:,:) :: dFHS
    COMPLEX(dp), ALLOCATABLE , DIMENSION(:,:,:,:) :: dFHS_k
    REAL(dp)   , ALLOCATABLE , DIMENSION(:,:,:,:) :: dFex ! gradient in real space
    CHARACTER(2) :: hs_functional ! hard sphere functional = PY for Percus-Yevick or CS for Carnahan-Starling
    REAL(dp), PARAMETER :: inv8pi  = 1.0_dp/(  8.0_dp * pi)
    REAL(dp), PARAMETER :: inv12pi = 1.0_dp/( 12.0_dp * pi)
    REAL(dp), PARAMETER :: inv18pi = 1.0_dp/( 18.0_dp * pi)
    REAL(dp), PARAMETER :: inv24pi = 1.0_dp/( 24.0_dp * pi)
    REAL(dp), PARAMETER :: inv36pi = 1.0_dp/( 36.0_dp * pi)

    IF (nb_species/=1) STOP "I stop because FMT calculations are only possible with one solvent species."

    CALL CPU_TIME ( time0 ) ! init timer

    nfft1 = spaceGrid%n_nodes(1)
    nfft2 = spaceGrid%n_nodes(2)
    nfft3 = spaceGrid%n_nodes(3)
    dV = spaceGrid%dv
    kT = thermocond%kbT

    allocate( rho (nfft1,nfft2,nfft3,nb_species) ,SOURCE=0._dp)
    icg = 0
    do s=1,size(solvent)
      do i=1,nfft1
        do j=1,nfft2
          do k=1,nfft3
            local_density=0._dp
            do o=1,angGrid%n_angles
              do p=1,molRotGrid%n_angles
                icg = icg + 1
                local_density = local_density +angGrid%weight(o)*cg_vect(icg)**2 *molRotGrid%weight(p)
              end do
            end do
            ! correct by *(8*pi²/n)**-1 as the integral over all orientations o and psi is 4pi and 2pi/n
            ! at the same time integrate rho in order to count the total number of implicit molecules. here we forget the integration factor = n_0 * dV
            rho(i,j,k,s) = local_density*molRotSymOrder/(fourpi*twopi)
          end do
        end do
      end do
    end do

    if ( all(rho<=epsilon(1._dp)) ) then
      Ffmt = -SUM( hs%Fexc0 * mole_fraction )
      FF=FF+Ffmt
      return
    end if

    ! total number of molecules of each species
    do concurrent ( s=1:nb_species )
        nb_molecules(s) = SUM( rho(:,:,:,s ) ) * solvent(s)%n0 * mole_fraction(s) * dV
    end do




    wdl=0 ! 0:3 for Kierliek Rosinberg
    wdu=3
    if( lbound(hs(1)%w_k,4) /= 0 ) then
      stop "problem with lbound l95 of energy_hard_sphere_fmt.f90"
    else if( ubound(hs(1)%w_k,4) /= 3 ) then
      stop "problem with ubound l98 of energy_hard_sphere_fmt.f90"
    else
      allocate( wd (nfft1,nfft2,nfft3, wdl:wdu) ,SOURCE=0._dp)
    end if

    do s=1,nb_species
      fftw3%in_forward = rho(:,:,:,s)
      call dfftw_execute ( fftw3%plan_forward ) ! fourier transform the density
      do i=wdl,wdu
        fftw3%in_backward = fftw3%out_forward * solvent(s)%n0 * mole_fraction(s) * hs(s)%w_k(:,:,:,i) ! rho_k(s) * w_k(i)
        call dfftw_execute( fftw3%plan_backward )
        wd(:,:,:,i) = wd(:,:,:,i) + fftw3%out_backward /(nfft1*nfft2*nfft3)
      end do
    end do
    deallocate(rho)



    ! check if the hard sphere functional is Percus-Yevick or Carnahan-Starling
    ! Get the free energy functional that should be used. For now Percus Yevick and Carnahan Starling only. May be expanded.
    hs_functional = input_char('hs_functional')
    if( hs_functional=='py') then
      hs_functional='PY'
    else if( hs_functional=='cs') then
      hs_functional='CS'
    end if

    if( hs_functional=='PY') then
      if( any(wd(:,:,:,3)>=1) ) then
        print*, "problem in log(1-w3)"
        call error_message_energy_hard_sphere_fmt (i,j,k,w0,w1,w2,w3)
        call process_output
      end if
      Ffmt =  -sum(wd(:,:,:,0)*LOG(1-wd(:,:,:,3))) &
            +sum(wd(:,:,:,1)*wd(:,:,:,2)/(1-wd(:,:,:,3))) &
            +inv24pi*sum(wd(:,:,:,2)**3/(1-wd(:,:,:,3))**2)

    else if( hs_functional=='CS') then
      if( any(wd(:,:,:,3)<=epsilon(1._dp)) ) then
        print*,"I stop before dividing by",w3
        call error_message_energy_hard_sphere_fmt (i,j,k,w0,w1,w2,w3)
        call process_output
      end if
      Ffmt = sum(wd(:,:,:,1)*wd(:,:,:,2)/(1-wd(:,:,:,3))) &
             +sum( (inv36pi*wd(:,:,:,2)**3 / wd(:,:,:,3)**2 -wd(:,:,:,0)) * LOG((1-wd(:,:,:,3)))  ) &
             +sum( inv36pi*wd(:,:,:,2)**3/(wd(:,:,:,3)*(1-wd(:,:,:,3))**2)  )

    end if

    Ffmt = Ffmt*kT*dV -SUM(hs%excchempot * nb_molecules) -SUM( hs%Fexc0 * mole_fraction )
    FF = FF + Ffmt

    ! gradients
    ! dFHS_i and weighted_density_j are arrays of dimension (nfft1,nfft2,nfft3)
    ! one_min_weighted_density_3 is dummy for speeding up and simplifying while using unnecessary memory.
    allocate ( one_min_wd_3 ( nfft1 , nfft2 , nfft3 ) ,SOURCE=1.0_dp-wd(:,:,:,3))

    ! excess free energy per atom derived wrt weighted density 0 (eq. A5 in Sears2003)
    allocate ( dFHS(nfft1,nfft2,nfft3,wdl:wdu) ,SOURCE=0._dp)

    ! Perkus Yevick
    IF ( hs_functional == 'PY' ) THEN
        IF( ANY(one_min_wd_3<=EPSILON(1.0_dp))) STOP "I found some value in one_min_wd_3 that cannot go in log"
        dFHS(:,:,:,0) = - log ( one_min_wd_3 )
        dFHS(:,:,:,1) = wd(:,:,:,2) / one_min_wd_3
        dFHS(:,:,:,2) = wd(:,:,:,1) / one_min_wd_3 + inv8pi * dFHS(:,:,:,1) ** 2
        dFHS(:,:,:,3) = ( wd(:,:,:,0) + wd(:,:,:,1) * dFHS(:,:,:,1) )/one_min_wd_3  + inv12pi*dFHS(:,:,:,1)**3

    ! Carnahan Starling
    ELSE IF ( hs_functional == 'CS' ) THEN

        IF( ANY(one_min_wd_3<=EPSILON(1.0_dp))) STOP "I found some value in one_min_wd_3 that cannot go in log"
        IF( ANY(ABS(wd(:,:,:,3))<=EPSILON(1.0_dp))) STOP "I found some value in w3 that would lead to divide by 0"

        dFHS(:,:,:,0) = - log ( one_min_wd_3 )

        dFHS(:,:,:,1) = wd(:,:,:,2) / one_min_wd_3

        dFHS(:,:,:,2) = - inv12pi * ( wd(:,:,:,2) / wd(:,:,:,3) ) ** 2 * dFHS(:,:,:,0) &
                + wd(:,:,:,1) / one_min_wd_3 + inv12pi * dFHS(:,:,:,1) ** 2 / wd(:,:,:,3)

        dFHS(:,:,:,3) = inv18pi * dFHS(:,:,:,0) * ( wd(:,:,:,2) / wd(:,:,:,3) ) ** 3 &
                - ( inv36pi * wd(:,:,:,2) ** 3 / wd(:,:,:,3) ** 2 - wd(:,:,:,0) ) /&
                                    one_min_wd_3 &
                + wd(:,:,:,1) * dFHS(:,:,:,1) / one_min_wd_3 + inv36pi * &
                wd(:,:,:,2) ** 3 / wd(:,:,:,3) ** 2 * &
                ( 3.0_dp * wd(:,:,:,3) - 1.0_dp ) / one_min_wd_3 ** 3
    END IF

    ! deallocate weighted_densities
    deallocate ( wd )
    deallocate ( one_min_wd_3 )

    ! compute gradients in k space
    allocate ( dFHS_k(nfft1/2+1,nfft2,nfft3,0:3) ,SOURCE=zeroC)

    ! FFT dFHS for computing convolution
    do i= wdl,wdu
        fftw3%in_forward = dFHS(:,:,:,i)
        CALL dfftw_execute ( fftw3%plan_forward )
        dFHS_k(:,:,:,i) = fftw3%out_forward
    end do

    ! deallocate useless
    deallocate ( dFHS )

    ! compute final gradient in k-space
    allocate ( dFex (nfft1,nfft2,nfft3,nb_species) ,SOURCE=0._dp)
    do s=1,nb_species
      fftw3%in_backward = SUM( dFHS_k * hs(s)%w_k ,4)
      CALL dfftw_execute ( fftw3%plan_backward )
      dFex(:,:,:,s) = fftw3%out_backward / (nfft1*nfft2*nfft3)
    end do
    deallocate ( dFHS_k )


    ! transfer in rank 1 vector dF
    icg = 0
    do s = 1 , nb_species
        do i = 1, nfft1
            do j = 1, nfft2
                do k = 1, nfft3
                    do o = 1, angGrid%n_angles
                        do p=1, molRotGrid%n_angles
                            icg = icg + 1
                            psi = cg_vect ( icg )
    ! AHAAAAAAAAAAAAATTENTION ICI LE 12 SEPTEMBRE 2011 JE ME RENDS COMPTE QU iL Y A PEUT ETRE UN FACTEUR WEIGHT(o) QUI MANQUE, CA DEPEND DE LA DEF DE N0 par RAPPORT à rho_0
                            dF ( icg ) = dF ( icg ) &
                            + 2.0_dp * psi * solvent(s)%rho0 * dV * ( kT * dFex(i,j,k,s) - hs(s)%excchempot )&
                            *angGrid%weight(o)*molRotGrid%weight(p)
    ! ATTENTION J'AI MULTIPLIE PAR WEIGHT(O) QD PASSAGE A angGrid%n_angles /=1 LE 18 JUILLET 2011
                        end do
                    end do
                end do
            end do
        end do
    end do
    deallocate (dFex)

    CALL CPU_TIME ( time1 )


    CONTAINS

    !===============================================================================================================================
    ! this SUBROUTINE prints error message related to SUBROUTINE excess_cs_hard_sphere
    ! it may stop program execution depending on the error.
    !===============================================================================================================================
    SUBROUTINE error_message_energy_hard_sphere_fmt ( i , j , k , w0 , w1 , w2 , w3 )
        USE precision_kinds ,ONLY: i2b, dp
        IMPLICIT NONE
        INTEGER(i2b), INTENT(IN) :: i,j,k
        REAL(dp), INTENT(IN) :: w0,w1,w2,w3
        PRINT*,'i , j , k = ',i,j,k
        PRINT*,'log (1-w3<=0) in energy_hard_sphere_fmt.f90. Critical stop'
        PRINT*,'w0, w1, w2, w3 = ' ,w0,w1,w2,w3
        STOP
    END SUBROUTINE error_message_energy_hard_sphere_fmt

END SUBROUTINE energy_hard_sphere_fmt
