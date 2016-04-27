!     Last change:  LB   14 Dec 2006    5:04 pm

module lecture
  !
  interface lis
    module PROCEDURE lis_entier,lis_entier_tableau, &
    lis_reel4,lis_reel4_tableau, &
    lis_reel8,lis_reel8_tableau
  end interface
  !
contains
  !
  subroutine lis_entier(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)
    INTEGER, INTENT(INOUT) :: i1
    INTEGER, INTENT(INOUT), optional:: i2,i3,i4,i5,i6,i7,i8,i9,i10
    INTEGER :: ient
    IF(.not.present(i2)) READ(*,*,ERR=100) ient
    IF(PRESENT(i2).and..not.present(i3)) READ(*,*,ERR=100) ient,i2
    IF(PRESENT(i3).and..not.present(i4)) READ(*,*,ERR=100) ient,i2,i3
    IF(PRESENT(i4).and..not.present(i5)) READ(*,*,ERR=100) ient,i2,i3,i4
    IF(PRESENT(i5).and..not.present(i6)) READ(*,*,ERR=100) ient,i2,i3,i4,i5
    IF(PRESENT(i6).and..not.present(i7)) READ(*,*,ERR=100) ient,i2,i3,i4,i5,i6
    IF(PRESENT(i7).and..not.present(i8)) READ(*,*,ERR=100) ient,i2,i3,i4,i5,i6,i7
    IF(PRESENT(i8).and..not.present(i9)) READ(*,*,ERR=100) ient,i2,i3,i4,i5,i6,i7,i8
    IF(PRESENT(i9).and..not.present(i10)) READ(*,*,ERR=100) ient,i2,i3,i4,i5,i6,i7,i8,i9
    IF(PRESENT(i10)) READ(*,*,ERR=100) ient,i2,i3,i4,i5,i6,i7,i8,i9,i10
    i1=ient
    return
    100 continue
    IF(.not.PRESENT(i2)) WRITE(*,*) i1
    IF(PRESENT(i2).and..not.present(i3)) WRITE(*,'(1x,i6)') i1,i2
    IF(PRESENT(i3).AND..not.present(i4)) WRITE(*,'(1x,i6)') i1,i2,i3
    IF(PRESENT(i4).and..not.present(i5)) WRITE(*,'(1x,i6)') i1,i2,i3,i4
    IF(PRESENT(i5).and..not.present(i6)) WRITE(*,'(1x,i6)') i1,i2,i3,i4,i5
    IF(PRESENT(i6).and..not.present(i7)) WRITE(*,'(1x,i6)') i1,i2,i3,i4,i5,i6
    IF(PRESENT(i7).and..not.present(i8)) WRITE(*,'(1x,i6)') i1,i2,i3,i4,i5,i6,i7
    IF(PRESENT(i8).and..not.present(i9)) WRITE(*,'(1x,i6)') i1,i2,i3,i4,i5,i6,i7,i8
    IF(PRESENT(i9).and..not.present(i10)) WRITE(*,'(1x,i6)') i1,i2,i3,i4,i5,i6,i7,i8,i9
    IF(PRESENT(i10)) WRITE(*,'(1x,i6)') i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
    WRITE(*,*)
  end subroutine
  !
  subroutine lis_entier_tableau(itab)
    INTEGER,INTENT(INOUT),DIMENSION(:) :: itab
    integer j, ient, n
    n=SIZE(itab)
    read(*,*,ERR=100) ient,(itab(j),j=2,n)
    itab(1)=ient
    return
    100 continue
    WRITE(*,'(1x,i6)') (itab(j),j=1,n)
    WRITE(*,*)
  end subroutine
  !
  subroutine lis_reel4(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
    REAL(4), INTENT(INOUT) :: x1
    REAL(4), INTENT(INOUT), optional :: x2,x3,x4,x5,x6,x7,x8,x9,x10
    REAL(4) :: reel
    IF(.not.present(x2)) READ(*,*,ERR=100) reel
    IF(PRESENT(x2).and..not.present(x3)) READ(*,*,ERR=100) reel,x2
    IF(PRESENT(x3).and..not.present(x4)) READ(*,*,ERR=100) reel,x2,x3
    IF(PRESENT(x4).and..not.present(x5)) READ(*,*,ERR=100) reel,x2,x3,x4
    IF(PRESENT(x5).and..not.present(x6)) READ(*,*,ERR=100) reel,x2,x3,x4,x5
    IF(PRESENT(x6).and..not.present(x7)) READ(*,*,ERR=100) reel,x2,x3,x4,x5,x6
    IF(PRESENT(x7).and..not.present(x8)) READ(*,*,ERR=100) reel,x2,x3,x4,x5,x6,x7
    IF(PRESENT(x8).and..not.present(x9)) READ(*,*,ERR=100) reel,x2,x3,x4,x5,x6,x7,x8
    IF(PRESENT(x9).and..not.present(x10)) READ(*,*,ERR=100) reel,x2,x3,x4,x5,x6,x7,x8,x9
    IF(PRESENT(x10)) READ(*,*,ERR=100) reel,x2,x3,x4,x5,x6,x7,x8,x9,x10
    x1=reel
    return
    100 continue
    IF(.not.PRESENT(x2)) WRITE(*,'(1x,e12.6)') x1
    IF(PRESENT(x2).and..not.present(x3)) WRITE(*,'(1x,e12.6)') x1,x2
    IF(PRESENT(x3).AND..not.present(x4)) WRITE(*,'(1x,e12.6)') x1,x2,x3
    IF(PRESENT(x4).and..not.present(x5)) WRITE(*,'(1x,e12.6)') x1,x2,x3,x4
    IF(PRESENT(x5).and..not.present(x6)) WRITE(*,'(1x,e12.6)') x1,x2,x3,x4,x5
    IF(PRESENT(x6).and..not.present(x7)) WRITE(*,'(1x,e12.6)') x1,x2,x3,x4,x5,x6
    IF(PRESENT(x7).and..not.present(x8)) WRITE(*,'(1x,e12.6)') x1,x2,x3,x4,x5,x6,x7
    IF(PRESENT(x8).and..not.present(x9)) WRITE(*,'(1x,e12.6)') x1,x2,x3,x4,x5,x6,x7,x8
    IF(PRESENT(x9).and..not.present(x10)) WRITE(*,'(1x,e12.6)') x1,x2,x3,x4,x5,x6,x7,x8,x9
    IF(PRESENT(x10)) WRITE(*,'(1x,e12.6)') x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
    WRITE(*,*)
  end subroutine
  !
  subroutine lis_reel4_tableau(xtab)
    REAL(4),INTENT(INOUT),DIMENSION(:) :: xtab
    real(4) reel
    integer n, j
    n=SIZE(xtab)
    read(*,*,ERR=100) reel,(xtab(j),j=2,n)
    xtab(1)=reel
    return
    100 continue
    WRITE(*,'(1x,e12.6)') (xtab(j),j=1,n)
    WRITE(*,*)
  end subroutine
  !
  subroutine lis_reel8(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
    REAL(8), INTENT(INOUT) :: x1
    REAL(8), INTENT(INOUT), optional :: x2,x3,x4,x5,x6,x7,x8,x9,x10
    REAL(8) :: reel
    IF(.not.present(x2)) READ(*,*,ERR=100) reel
    IF(PRESENT(x2).and..not.present(x3)) READ(*,*,ERR=100) reel,x2
    IF(PRESENT(x3).and..not.present(x4)) READ(*,*,ERR=100) reel,x2,x3
    IF(PRESENT(x4).and..not.present(x5)) READ(*,*,ERR=100) reel,x2,x3,x4
    IF(PRESENT(x5).and..not.present(x6)) READ(*,*,ERR=100) reel,x2,x3,x4,x5
    IF(PRESENT(x6).and..not.present(x7)) READ(*,*,ERR=100) reel,x2,x3,x4,x5,x6
    IF(PRESENT(x7).and..not.present(x8)) READ(*,*,ERR=100) reel,x2,x3,x4,x5,x6,x7
    IF(PRESENT(x8).and..not.present(x9)) READ(*,*,ERR=100) reel,x2,x3,x4,x5,x6,x7,x8
    IF(PRESENT(x9).and..not.present(x10)) READ(*,*,ERR=100) reel,x2,x3,x4,x5,x6,x7,x8,x9
    IF(PRESENT(x10)) READ(*,*,ERR=100) reel,x2,x3,x4,x5,x6,x7,x8,x9,x10
    x1=reel
    return
    100 continue
    IF(.not.PRESENT(x2)) WRITE(*,'(1x,d23.15)') x1
    IF(PRESENT(x2).and..not.present(x3)) WRITE(*,'(1x,d23.15)') x1,x2
    IF(PRESENT(x3).AND..not.present(x4)) WRITE(*,'(1x,d23.15)') x1,x2,x3
    IF(PRESENT(x4).and..not.present(x5)) WRITE(*,'(1x,d23.15)') x1,x2,x3,x4
    IF(PRESENT(x5).and..not.present(x6)) WRITE(*,'(1x,d23.15)') x1,x2,x3,x4,x5
    IF(PRESENT(x6).and..not.present(x7)) WRITE(*,'(1x,d23.15)') x1,x2,x3,x4,x5,x6
    IF(PRESENT(x7).and..not.present(x8)) WRITE(*,'(1x,d23.15)') x1,x2,x3,x4,x5,x6,x7
    IF(PRESENT(x8).and..not.present(x9)) WRITE(*,'(1x,d23.15)') x1,x2,x3,x4,x5,x6,x7,x8
    IF(PRESENT(x9).and..not.present(x10)) WRITE(*,'(1x,d23.15)') x1,x2,x3,x4,x5,x6,x7,x8,x9
    IF(PRESENT(x10)) WRITE(*,'(1x,d23.15)') x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
    WRITE(*,*)
  end subroutine
  !
  subroutine lis_reel8_tableau(xtab)
    REAL(8),INTENT(INOUT),DIMENSION(:) :: xtab
    REAL(8) :: reel
    integer j, n
    n=SIZE(xtab)
    read(*,*,ERR=100) reel,(xtab(j),j=2,n)
    xtab(1)=reel
    return
    100 continue
    WRITE(*,'(1x,d23.15)') (xtab(j),j=1,n)
    WRITE(*,*)
  end subroutine
  !
end module
