program test
    use iso_c_binding, only: dp=>C_DOUBLE
    implicit none
    integer :: mmax
    integer :: i, alp00000, j
    integer, parameter :: nalp(5)=[6,75,252,877,2002]
    integer, dimension(:), allocatable :: alp, m, n, mu, nu, khi
    character(3) :: tmpchar
    real(dp) :: k
    complex(dp), dimension(:), allocatable :: cmnmunukhi
    do mmax=1,5
        j=nalp(mmax)
        allocate(alp(j), m(j), n(j), mu(j), nu(j), khi(j) )
        allocate(cmnmunukhi(j))
        select case (mmax)
        case(1)
            open(10,file="ck_nonzero_nmax1_ml")
            open(11,file="cs_from_luc_mmax1")
        case(2)
            open(10,file="ck_nonzero_nmax2_ml")
            open(11,file="cs_from_luc_mmax2")
        case(3)
            open(10,file="ck_nonzero_nmax3_ml")
            open(11,file="cs_from_luc_mmax3")
        case(4)
            open(10,file="ck_nonzero_nmax4_ml")
            open(11,file="cs_from_luc_mmax4")
        case(5)
            open(10,file="ck_nonzero_nmax5_ml")
            open(11,file="cs_from_luc_mmax5")
        end select
        do i=1,10
            read(10,*)
        end do
        read(10,*) tmpchar, alp
        read(10,*) tmpchar, m
        read(10,*) tmpchar, n
        read(10,*) tmpchar, mu
        read(10,*) tmpchar, nu
        read(10,*) tmpchar, khi
        read(10,*) tmpchar
        do i=1,j
            if (m(i)==0.and.n(i)==0.and.mu(i)==0.and.nu(i)==0.and.khi(i)==0) then
                alp00000=i
            end if
        end do
        do i=1,199
            read(10,*) k, cmnmunukhi
            write(11,*) k, real(cmnmunukhi(alp00000))
        end do
        close(10)
        close(11)
        deallocate(alp,m,n,mu,nu,khi)
        deallocate(cmnmunukhi)
    end do
end program test
