! This SUBROUTINE computes the radial distribution function around each solute site
SUBROUTINE compute_rdf ( array , filename )

    USE precision_kinds, ONLY: dp , i2b
    USE system, ONLY: x_mol , y_mol , z_mol , id_mol , id_solv, nb_species, spaceGrid, soluteSite
    USE input, ONLY: input_dp, input_int
    
    IMPLICIT NONE

    REAL(dp), DIMENSION(spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3),nb_species), INTENT(IN) :: array
    CHARACTER( 50 ) , intent(in) :: filename
    INTEGER(i2b):: nb_id_solv ! Total number of solvent site kinds
    INTEGER(i2b):: nb_id_mol ! Total number of solute site kinds
    INTEGER(i2b):: nbins !> Total number of bins in which rdf is histogrammed
    REAL(dp):: RdfMaxRange !> Maximum distance to which calculate radial distribution function
    REAL(dp):: delta_r !> Width of each bin
    REAL(dp):: r !> Distance ^2 between each of the n solute sites and each grid point i,j,k
    REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: rdf
    INTEGER(i2b), ALLOCATABLE, DIMENSION (:,:) :: recurrence_bin
    INTEGER(i2b):: i , j , k , n , bin, s
    REAL(dp):: xnm2 , ynm2 , znm2 ! dummy

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
delta_r = RdfMaxRange / REAL ( nbins , dp ) ! Width of each bin of the histogram
! open the file
open ( 10 , file = filename , form = 'formatted' )
100 format ( 2 ( xF10.5 ) )
write ( 10 , * ) '# r  rdf'
! prepare histogram
allocate ( recurrence_bin ( SIZE(soluteSite), 0:nbins ) )
allocate ( rdf ( SIZE(soluteSite), 0:nbins ) )

! FOR EACH SPECIES
do s = 1 , nb_species
  recurrence_bin = 0
  rdf = 0.0_dp
  write ( 10 , * ) '# species number ', s
  !> Transform array(position) in rdf(radialdistance)
  ! counts the total number of appearence of a value in each bin
  do n = 1 , SIZE(soluteSite)
    do i = 1 , spaceGrid%n_nodes(1)
      xnm2 = ( REAL ( i - 1 , dp ) * spaceGrid%dl(1) - x_mol ( n ) ) ** 2
      do j = 1 , spaceGrid%n_nodes(2)
        ynm2 = ( REAL ( j - 1 , dp ) * spaceGrid%dl(2) - y_mol ( n ) ) ** 2
        do k = 1 , spaceGrid%n_nodes(3)
          znm2 = ( REAL ( k - 1 , dp ) * spaceGrid%dl(3) - z_mol ( n ) ) ** 2
          r = sqrt ( xnm2 + ynm2 + znm2 ) ! Distance r between each solute site N to each grid point I,J,K
          if ( r > RdfMaxRange ) cycle !** 2 ) cycle
          bin = int ( r / delta_r + 0.5_dp ) !> In which bin of histogram are we.
          if ( bin > nbins ) cycle ! out of range
          recurrence_bin ( n , bin ) = recurrence_bin ( n , bin ) + 1
          rdf ( n , bin ) = rdf ( n , bin ) + array ( i , j , k , s )
        END DO
      END DO
    END DO
  END DO

  !> Average rdf on each bin
  where ( recurrence_bin /= 0 )
     rdf = rdf / REAL ( recurrence_bin , dp )
  ELSEWHERE
     rdf = 0.0_dp
  end where

  ! Write rdf in output folder
  do n = 1 , SIZE(soluteSite)
     write ( 10 , * ) '#site ' , n
     do bin = 0 , nbins
           write ( 10 , 100 ) (REAL ( bin , dp )) * delta_r , rdf ( n , bin ) 
     END DO
     write ( 10 , * )
  END DO
  ! white line for xmgrace
  write ( 10 , * )
END DO
close(10)

deallocate ( recurrence_bin )
deallocate ( rdf )

END SUBROUTINE compute_rdf
