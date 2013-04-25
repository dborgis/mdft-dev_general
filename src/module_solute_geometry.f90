module solute_geometry
    use precision_kinds
    use system, only: nb_solute_sites, x_mol, y_mol, z_mol ! TODO should be generalized. nb_solute_sites is useless as it is the dimension of x_mol. isLinear and isPlanar could be usefull for any coordinates, not only the ones of the solute. the coordinates should thus be intent(in) of the functions.
    implicit none
    private
    public :: isLinear, isPlanar
    contains
    
        ! Test if the solute is linear
        logical pure function isLinear()
            integer:: n,sim1,sim2,sim3 !> dummy
            ! One site cannot be linear, while two are always, else ...
            if (nb_solute_sites==1) then
                islinear = .false.
            else if (nb_solute_sites==2) then
                islinear = .true.
            else
                ! Are there two common coordinates to all sites ?
                sim1=0
                sim2=0
                sim3=0
                do n=1,nb_solute_sites
                    if (x_mol(n)==y_mol(n)) sim1=sim1+1
                    if (x_mol(n)==z_mol(n)) sim2=sim2+1
                    if (y_mol(n)==z_mol(n)) sim3=sim3+1
                end do
                if (sim1==nb_solute_sites .or. sim2==nb_solute_sites .or. sim3==nb_solute_sites) then
                    islinear = .true.
                else
                    islinear = .false.
                end if
            end if
        end function
        
        
        ! Check solute planarity
        logical pure function isPlanar()
            real(dp), dimension(2:nb_solute_sites):: xvec,yvec,zvec !> vector between solute sites 1 and N (2<N<nb_solute_sites)
            integer :: n !> dummy
            !> Can be planar only if number of solute sites is 3 or more
            if (nb_solute_sites==1 .or. nb_solute_sites==2) then
                isPlanar = .false.
            else
                !> Compute all vector coordinates between first site and every others
                forall (n=2:nb_solute_sites)
                    xvec(n) = x_mol(1)-x_mol(n)
                    yvec(n) = y_mol(1)-y_mol(n)
                    zvec(n) = z_mol(1)-z_mol(n)
                end forall
                !> Check if one coordinate is always the same (which is the case if we're in a plan)
                if ( minval(xvec)==maxval(xvec) .or. minval(yvec)==maxval(yvec) .or. minval(zvec)==maxval(zvec)) then
                    isPlanar = .true.
                else
                    isPlanar = .false.
                end if
            end if
        end function

end module
