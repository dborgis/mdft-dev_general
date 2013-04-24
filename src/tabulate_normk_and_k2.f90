! This subroutine tabulates the value of norm(k) for each vector k (l,m,n) ((nfft1,nfft2,nfft3)) and its square

subroutine tabulate_normk_and_k2



use precision_kinds , only : i2b , dp

use system , only : nfft1 , nfft2 , nfft3 , Lx , Ly , Lz

use fft , only : norm_k , k2 , kx , ky , kz , k2_nocoef! what we want to tabulate

use constants , only : twopi



implicit none



integer ( kind = i2b ) :: l , m , n , m1 , m2 , m3 , nf1 , nf2 , nf3 ! dummy

real ( kind = dp ) :: kx2 , ky2 , kz2 ! norm squared of projection of k in each direction

real ( kind = dp ) :: twopioverlx , twopioverly , twopioverLz





allocate ( norm_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )

allocate ( k2 ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )

allocate ( k2_nocoef ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )

allocate ( kx ( nfft1 / 2 + 1 ) )

allocate ( ky ( nfft2 ) )

allocate ( kz ( nfft3 ) )

! init dummy

nf1 = nfft1 / 2

nf2 = nfft2 / 2

nf3 = nfft3 / 2

twopioverlx = twopi / Lx

twopioverly = twopi / Ly

twopioverLz = twopi / Lz


! loop over nfft1, 2 & 3 paying attention to iner loop being 1 and outer loop 3

do n = 1 , nfft3

  if ( n <= nf3 ) then

    m3 = n - 1

  else

    m3 = n - 1 - nfft3

  end if

  kz ( n ) = twopioverlz * real ( m3 , kind = dp )

  kz2 = kz ( n ) ** 2

  do m = 1 , nfft2

    if ( m <= nf2 ) then

      m2 = m - 1

    else

      m2 = m - 1 - nfft2

    end if

    ky ( m ) = twopioverly * real ( m2 , kind = dp )

    ky2 = ky ( m ) ** 2

    do l = 1 , nf1 + 1

      if ( l <= nf1 ) then

        m1 = l - 1

      else

        m1 = l - 1 - nfft1

      end if

      kx ( l ) = twopioverlx * real ( m1 , kind = dp )

      kx2 = kx ( l ) ** 2



      norm_k ( l , m , n ) = sqrt ( kx2 + ky2 + kz2 )

      k2 ( l , m , n ) = kx2 + ky2 + kz2

      k2_nocoef ( l , m , n ) = real( m1 ,kind=dp) ** 2 + real( m2 ,kind=dp) ** 2 + real( m3 ,kind=dp) ** 2

    end do

  end do

end do


end subroutine tabulate_normk_and_k2
