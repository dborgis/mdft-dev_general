! this subroutine gets the informations in input file and then allocate, prepare, compute
! needed data

subroutine prepare_minimizer


use precision_kinds , only : i2b , dp

use input , only : input_line, input_int, input_dp,input_char

use cg , only : nbd , iwa , ll , uu , wa , dF , cg_vect , epsmch , factr , epsg , ncg , mcg , iprint , minimizer_type ,&
                itermax , minimizer_iter , pgtol

use system , only : nfft1 , nfft2 , nfft3 , nb_omega , nb_species , nb_psi

implicit none

integer ( kind = i2b ) :: i , j ! dummy



! total number of variables to optimize

ncg = nfft1 * nfft2 * nfft3 * nb_omega *nb_psi* nb_species

! get maximum number of iterations itermax
itermax=input_int('maximum_iteration_nbr')

write (*,*) 'itermax = ' , itermax

! get epsg
epsg=input_dp('epsg')
write (*,*) 'epsg = ' , epsg


! get pgtol
pgtol=input_dp('pgtol')
write (*,*) 'pgtol = ' , pgtol

! get from input the type of minimization one wants to do

minimizer_type=trim(adjustl(input_char('minimizer')))

write (*,*) 'minimizer is ' , minimizer_type (1:4)

! init minimizer iteration counter

minimizer_iter = 0

! init X

allocate ( cg_vect ( ncg ) )

! init grad(F)

allocate ( dF ( ncg ) )

! if minimizer is bfgs then init what has to be init

if( minimizer_type(1:4) == 'bfgs') then

  mcg = 4 ! dont understand EXPLIQUER

  allocate ( ll(ncg) ) ! lower bound of cg_vect

  ll = 0.0_dp ! init

  allocate ( uu(ncg) ) ! uppder bound of cg_vect

  uu = 10.0_dp ! init

  allocate ( nbd(ncg) )

  nbd = 0 ! nbd(i) is : 0 for no constaint, 1 for lower bounded by ll, 3 for upper bounded by uu, 2 for lower and upper bounded

  epsmch = epsilon(1.0_dp)  !  Precision de la machine

  factr = epsg / epsmch ! convergence criteria over energy

  iprint = 1

  allocate ( iwa ( 3 * ncg ) )

  allocate ( wa ( 2 * mcg * ncg + 4 * ncg + 11 * mcg **2 + 8 * mcg ) )

end if


 

end subroutine prepare_minimizer
