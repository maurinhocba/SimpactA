      SUBROUTINE SORTND(A,N)
      IMPLICIT NONE
      INTEGER (kind=4) :: A(*),N

!***  SUPPLIED BY DON CALAHAN

!   PARTITION SORTING ALGORITHM
! REFERENCE COLLECTED ALGORITHMS OF THE ACM - 63,64,65

      INTEGER (kind=4) :: IHIGH(32),ILOW(32),NSEGS,IL,IH,ISEPX,ISEP,IXL, &
     &                    IXH,IT
! INITIALIZE
      NSEGS = 1
      IL = 1
      IH = N
! IF NO ELEMENTS IN THIS SEGMENT DO NOTHING
   10 IF (IL >= IH) GO TO 80
! CHOOSE ISEP (SEPARATION ENTRY):
!  MAKE A(IL) <= A((IL+IH)/2) <= A(IH) BY INTERCHANGE
!  SET ISEP= A((IL+IH)/2)
   20 ISEPX = (IH + IL) / 2
      ISEP = A(ISEPX)
! IXL IS LOWER SEGMENT INDEX (CURRENT)
      IXL = IL
! MAKE A(IL) <= A(ISEPX)
      IF (A(IL) <= ISEP) GO TO 30
      A(ISEPX) = A(IL)
      A(IL) = ISEP
      ISEP = A(ISEPX)
! IXH IS HIGHEST SEGMENT INDEX (CURRENT)
   30 IXH = IH
! MAKE A(IH) >= A(ISEPX)
      IF (A(IH) >= ISEP) GO TO 50
      A(ISEPX) = A(IH)
      A(IH) = ISEP
      ISEP = A(ISEPX)
! MAKE A(IL) <= A(ISEPX)
      IF (A(IL) <= ISEP) GO TO 50
      A(ISEPX) = A(IL)
      A(IL) = ISEP
      ISEP = A(ISEPX)
      GO TO 50
! EXCHANGE LOW PART ENTRY WHICH IS GREATER THAN SEPARATOR WITH HIGH
! PART ENTRY WHICH IS LESS THAN OR EQUAL TO THE SEPARATOR VALUE.
   40 IT = A(IXH)
      A(IXH) = A(IXL)
      A(IXL) = IT
! MOVE DOWN UPPER SEGMENT AS FAR AS WE CAN
   50 IXH = IXH - 1
      IF (A(IXH) > ISEP) GO TO 50
! MOVE UP LOWER SEGMENT AS FAR AS WE CAN
   60 IXL = IXL + 1
      IF (A(IXL) < ISEP) GO TO 60
! NOTHING TO DO IF BOTH SEGMENTS HAVE AT MOST ONE ENTRY IN COMMON
      IF (IXL <= IXH) GO TO 40
! IF BOTH SEGMENTS OVERLAP THEN THEY ARE SEPARATED
! IN THIS CASE CONTINUE WITH SHORTER SEGMENT, STORING THE LONGER
      IF (IXH - IL <= IH - IXL) GO TO 70
! LOWER SEGMENT LONGER, CONTIN WITH UPPER AFTER SAVING LOWER
      ILOW(NSEGS) = IL
      IHIGH(NSEGS) = IXH
      IL = IXL
      NSEGS = NSEGS + 1
      GO TO 90
! UPPER SEGMENT LONGER, CONTIN WITH LOWER AFTER SAVING UPPER
   70 ILOW(NSEGS) = IXL
      IHIGH(NSEGS) = IH
      IH = IXH
      NSEGS = NSEGS + 1
      GO TO 90
! GET ANOTHER SEGMENT FOR PROCESSING IF THERE ARE ANY MORE
   80 NSEGS = NSEGS - 1
      IF (NSEGS == 0) RETURN
      IL = ILOW(NSEGS)
      IH = IHIGH(NSEGS)
! CONTINUE TO SEGMENT AS LONG AS LENGTH IS GREATER THAN 11
   90 IF (IH - IL >= 11) GO TO 20
      IF (IL == 1) GO TO 10
      GO TO 110
! SORT ELEMENTS WITHIN SEGMENT BY INTERCHANGE OF ADJACENT PAIRS
  100 IL = IL + 1
  110 IF (IL == IH) GO TO 80
      ISEP = A(IL + 1)
      IF (A(IL) <= ISEP) GO TO 100
      IXL = IL
  120 A(IXL + 1) = A(IXL)
      IXL = IXL - 1
      IF (ISEP < A(IXL)) GO TO 120
      A(IXL + 1) = ISEP
      GO TO 100
      END
!      SPAGHETTI FORTRAN, IMPOSSIBLE TO UNDERSTAND
!      SUBROUTINE sortnd(a,n)
!      !***  supplied by Don Calahan
!
!      !   partition sorting algorithm
!      ! reference Collected Algorithms of the ACM - 63,64,65
!      IMPLICIT NONE
!      ! dummy arguments
!      INTEGER (kind=4) :: a(*),n
!      ! local varibles
!      INTEGER (kind=4) :: ihigh(32),ilow(32),nsegs,il,ih,isepx,isep,ixl, &
!                          ixh,it,il1
!      ! initialize
!      nsegs = 1
!      il = 1
!      ih = n
!      ! IF no elements in this segment do nothing
!   10 IF (il >= ih) GO TO 80
!      ! choose isep (separation entry):
!      !  make a(il) <= a((il+ih)/2) <= a(ih) by interchange
!      !  set isep= a((il+ih)/2)
!   20 isepx = (ih + il) / 2
!      isep = a(isepx)
!      ! ixl is lower segment index (current)
!      ixl = il
!      ! make a(il) <= a(isepx)
!      IF (a(il) > isep) THEN
!        a(isepx) = a(il)
!        a(il) = isep
!        isep = a(isepx)
!      END IF
!      ! ixh is highest segment index (current)
!      ixh = ih
!      ! make a(ih) >= a(isepx)
!      IF (a(ih) < isep) THEN
!        a(isepx) = a(ih)
!        a(ih) = isep
!        isep = a(isepx)
!      END IF
!      ! make a(il) <= a(isepx)
!      IF (a(il) > isep) THEN
!        a(isepx) = a(il)
!        a(il) = isep
!        isep = a(isepx)
!      END IF
!      ! exchange low part entry which is greater than separator with high
!      ! part entry which is less than or equal to the separator value.
!      DO
!        ! move down upper segment as far as we can
!        DO
!          ixh = ixh - 1
!          IF (a(ixh) <= isep) EXIT
!        END DO
!        ! move up lower segment as far as we can
!        DO
!          ixl = ixl + 1
!          IF (a(ixl) >= isep) EXIT
!        END DO
!        ! nothing to do IF both segments have at most one entry in common
!        IF (ixl > ixh) EXIT
!        it = a(ixh)
!        a(ixh) = a(ixl)
!        a(ixl) = it
!      END DO
!      ! IF both segments overlap then they are separated
!      ! in this case continue with shorter segment, storing the longer
!      IF (ixh - il > ih - ixl) THEN
!        ! upper segment longer, contin with lower after saving upper
!        ilow(nsegs) = ixl
!        ihigh(nsegs) = ih
!        ih = ixh
!        nsegs = nsegs + 1
!      ELSE
!        ! lower segment longer, contin with upper after saving lower
!        ilow(nsegs) = il
!        ihigh(nsegs) = ixh
!        il = ixl
!        nsegs = nsegs + 1
!      END IF
!      ! continue to segment as long as length is greater than 11
!      DO
!        IF (ih - il >= 11) GO TO 20
!        IF (il == 1) GO TO 10
!        ! sort elements within segment by interchange of adjacent pairs
!        DO
!          IF (il == ih) EXIT
!          isep = a(il1)
!          IF (a(il) <= isep) CYCLE
!          ixl = il
!          DO
!            a(ixl + 1) = a(ixl)
!            ixl = ixl - 1
!            IF (isep >= a(ixl)) EXIT
!          END DO
!          a(ixl + 1) = isep
!          il = il1
!        END DO
!        ! get another segment for processing IF there are any more
!   80   nsegs = nsegs - 1
!        IF (nsegs == 0) EXIT
!        il = ilow(nsegs)
!        ih = ihigh(nsegs)
!      END DO
!
!      RETURN
!      END SUBROUTINE sortnd
