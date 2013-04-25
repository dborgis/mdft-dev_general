! This subroutine computes the radial distribution function around each solute site
subroutine compute_rdf ( array , filename )
use precision_kinds , only : dp , i2b
use system , only : nfft1 , nfft2 , nfft3 , deltax , deltay , deltaz , nb_solute_sites , x_mol , y_mol , z_mol , id_mol , id_solv, &
                    nb_species
use input , only : input_line, input_dp, input_int
implicit none
real(dp), intent(in) , dimension ( nfft1 , nfft2 , nfft3 , nb_species ) :: array
character ( 50 ) , intent(in) :: filename
integer(i2b):: nb_id_solv ! Total number of solvent site kinds
integer(i2b):: nb_id_mol ! Total number of solute site kinds
integer(i2b):: nbins !> Total number of bins in which rdf is histogrammed
real(dp):: RdfMaxRange !> Maximum distance to which calculate radial distribution function
real(dp):: delta_r !> Width of each bin
real(dp):: r !> Distance ^2 between each of the n solute sites and each grid point i,j,k
real(dp), allocatable , dimension ( : , : ) :: rdf
integer(i2b), allocatable , dimension ( : , : ) :: recurrence_bin
integer(i2b):: i , j , k , n , bin
real(dp):: xnm2 , ynm2 , znm2 ! dummy
integer(i2b):: species ! dummy between 1 and nb_species
!> Get the total number of solute kinds and solvent kinds
!TODO pour l'instant ne sert a rien mais ensuite ce n'est plus chaque site de solute qui doit compter mais chaque type de site
! on devra faire une moyenne sur tous les sites d'un meme type
! c'est pas bien compliqué à faire
! faut rajouter au bout de la routine, avant le write dans le fichier g.rdf, un test sur chaque solute
! si c'est le meme qu'un autre faut faire la moyenne et c'est tout
nb_id_solv = maxval ( id_solv )
nb_id_mol = maxval ( id_mol )
! Read the number of histogram steps
nbins=input_int('nbinsrdf')
! Read the maximum range of the radial distribution functions
RdfMaxRange=input_dp('rdfmaxrange')
!vRdfMaxRange=min(Lx,Ly,Lz)*sqrt(3.0_dp) !> Maximum distance to which calculate radial distribution function
delta_r = RdfMaxRange / real ( nbins , dp ) ! Width of each bin of the histogram
! open the file
open ( 10 , file = filename , form = 'formatted' )
100 format ( 2 ( xF10.5 ) )
write ( 10 , * ) '# r  rdf'
! prepare histogram
allocate ( recurrence_bin ( nb_solute_sites , 0 : nbins ) )
allocate ( rdf ( nb_solute_sites , 0 : nbins ) )
! FOR EACH SPECIES
do species = 1 , nb_species
  recurrence_bin = 0
  rdf = 0.0_dp
  write ( 10 , * ) '# species number ', species
  !> Transform array(position) in rdf(radialdistance)
  ! counts the total number of appearence of a value in each bin
  do n = 1 , nb_solute_sites
    do i = 1 , nfft1
      xnm2 = ( real ( i - 1 , dp ) * deltax - x_mol ( n ) ) ** 2
      do j = 1 , nfft2
        ynm2 = ( real ( j - 1 , dp ) * deltay - y_mol ( n ) ) ** 2
        do k = 1 , nfft3
          znm2 = ( real ( k - 1 , dp ) * deltaz - z_mol ( n ) ) ** 2
          r = sqrt ( xnm2 + ynm2 + znm2 ) ! Distance r between each solute site N to each grid point I,J,K
          if ( r > RdfMaxRange ) cycle !** 2 ) cycle
          bin = int ( r / delta_r + 0.5_dp ) !> In which bin of histogram are we.
          if ( bin > nbins ) cycle ! out of range
          recurrence_bin ( n , bin ) = recurrence_bin ( n , bin ) + 1
          rdf ( n , bin ) = rdf ( n , bin ) + array ( i , j , k , species )
        end do
      end do
    end do
  end do
  !> Average rdf on each bin
  where ( recurrence_bin /= 0 )
     rdf = rdf / real ( recurrence_bin , dp )
  elsewhere
     rdf = 0.0_dp
  end where
  ! Write rdf in output folder
  do n = 1 , nb_solute_sites
     write ( 10 , * ) '#site ' , n
     do bin = 0 , nbins
           write ( 10 , 100 ) (real ( bin , dp )) * delta_r , rdf ( n , bin ) 
     end do
     write ( 10 , * )
  end do
  ! white line for xmgrace
  write ( 10 , * )
end do
deallocate ( recurrence_bin )
deallocate ( rdf )
close(10)
write(*,*)filename,' written'
end subroutine compute_rdf
