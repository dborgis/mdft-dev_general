module solute_geometry
    USE precision_kinds
    use system, only: nb_solute_sites, solute ! TODO should be generalized. nb_solute_sites is useless as it is the dimension of x_mol. isLinear and isPlanar could be usefull for any coordinates, not only the ones of the solute. the coordinates should thus be intent(in) of the functions.
    IMPLICIT NONE
    private
    public :: isLinear, isPlanar
    contains
    
        ! Test if the solute is linear
        logical pure function isLinear()
            integer:: n,sim1,sim2,sim3 !> dummy
            ! One site cannot be linear, while two are always, ELSE ...
            if (nb_solute_sites==1) then
                islinear = .false.
            ELSE IF (nb_solute_sites==2) then
                islinear = .true.
            ELSE
                ! Are there two common coordinates to all sites ?
                sim1=0
                sim2=0
                sim3=0
                do n=1,nb_solute_sites
                    if (solute%site(n)%r(1)==solute%site(n)%r(2)) sim1=sim1+1
                    if (solute%site(n)%r(1)==solute%site(n)%r(3)) sim2=sim2+1
                    if (solute%site(n)%r(2)==solute%site(n)%r(3)) sim3=sim3+1
                END DO
                if (sim1==nb_solute_sites .or. sim2==nb_solute_sites .or. sim3==nb_solute_sites) then
                    islinear = .true.
                ELSE
                    islinear = .false.
                END IF
            END IF
        end function
        
        
        ! Check if solute has planar symetry
        logical pure function isPlanar()
            real(dp), dimension(2:nb_solute_sites):: xvec,yvec,zvec !> vector between solute sites 1 and N (2<N<nb_solute_sites)
            integer :: n !> dummy
            !> Can be planar only if number of solute sites is 3 or more
            if (nb_solute_sites==1 .or. nb_solute_sites==2) then
                isPlanar = .false.
            ELSE
                !> Compute all vector coordinates between first site and every others
                forall (n=2:nb_solute_sites)
                    xvec(n) = solute%site(1)%r(1)-solute%site(n)%r(1)
                    yvec(n) = solute%site(1)%r(2)-solute%site(n)%r(2)
                    zvec(n) = solute%site(1)%r(3)-solute%site(n)%r(3)
                end forall
                !> Check if one coordinate is always the same (which is the case if we're in a plan)
                if ( minval(xvec)==maxval(xvec) .or. minval(yvec)==maxval(yvec) .or. minval(zvec)==maxval(zvec)) then
                    isPlanar = .true.
                ELSE
                    isPlanar = .false.
                END IF
            END IF
        end function

end module
