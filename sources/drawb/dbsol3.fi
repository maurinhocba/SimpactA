SUBROUTINE dbsol3 (xl,      nrtm, ldbs,     irect,  &
                    fmax,   elmodul, nstatnew, nstatold, &
                    ndbpoin,  fdat,daz,  lmsr,    msr, nseg, &
                    nmn, ihflag,tdbf,nnode, fcont, x)
 !-----------------------------------------------------------------
 !  selects the sheet elements for the verification of intersection
 !-----------------------------------------------------------------

 IMPLICIT NONE

 !   d u m m y   a r g u m e n t s

 INTEGER, INTENT (IN) :: nrtm
 INTEGER :: ndbpoin
 INTEGER, INTENT (IN) :: nmn
 INTEGER :: ihflag
 INTEGER :: nnode
 REAL (kind=8) :: fmax
 REAL (kind=8) :: elmodul
 INTEGER :: ldbs (:, :)
 INTEGER :: irect (:, :)
 INTEGER :: nstatnew (nrtm)
 INTEGER :: nstatold (nrtm)
 INTEGER, INTENT (IN) :: lmsr (:)
 INTEGER, INTENT (IN) :: msr (:)
 INTEGER, INTENT (IN) :: nseg (:)
 REAL (kind=8) :: fcont (:, :)
 REAL (kind=8) :: x (:,:)
 REAL (kind=8) :: xl (:, :)
 REAL (kind=8) :: fdat (:, :)
 REAL (kind=8) :: tdbf (4)
 LOGICAL,INTENT(INOUT):: daz(:)

 !-----------------------------------------------
 !   l o c a l   v a r i a b l e s

 INTEGER :: i, ii, ind, jj, ll, num, npt, lpp
 REAL (kind=8) :: foldh
 !-----------------------------------------------

 nstatnew = - 1
 DO i = 1, nrtm
    IF (nstatold(i) == 1) THEN
       foldh = fdat (3, i)

       !  previously intersected elements are checked

!#if UNIX
!            CALL chquad (x, xl, nstatnew(i), nstatold(i), ldbs(:, i),   &
!     &                   irect(:, i), ndbpoin, fdat(:,i), daz(i),elmodul, fmax, &
!     &                   ndofn, fcont, tdbf, ihflag, foldh, &
!     &                   nnode, npoin)
!#else
        CALL chquad (x, xl, nstatnew(i), nstatold(i), ldbs(1:2, i),   &
                     irect(1:4,i),ndbpoin,fdat(1:3,i),daz(i),elmodul,fmax, &
                     fcont, tdbf, ihflag, foldh,  nnode)
!#endif
        DO ii = 1, 4
           ind = irect (ii, i)
           DO jj = 1, nmn
              IF (msr(jj) == ind) THEN
                 num = nseg (jj+1) - nseg (jj)
                 npt = nseg (jj) - 1
                 DO lpp = 1, num
                    ll = lmsr (lpp+npt)

                        ! and their neighbours
!#if UNIX
!                        IF (nstatold(ll) == 0 .AND. nstatnew(ll) ==(-1))
!     &                    CALL chquad
!     &                     (x, xl, nstatnew(ll),nstatold(ll),ldbs(:,ll),
!     &                      irect(:,ll),ndbpoin,fdat(:,ll),daz(i),elmodul,fmax,
!     &                      ndofn, fcont, tdbf, ihflag, foldh,
!     &                      nnode, npoin)
!#else
                    IF (nstatold(ll) == 0 .AND. nstatnew(ll) ==(-1)) &
                      CALL chquad &
                       (x, xl, nstatnew(ll),nstatold(ll),ldbs(:,ll), &
                        irect(:,ll),ndbpoin,fdat(:,ll),daz(i),elmodul,fmax, &
                        fcont, tdbf, ihflag, foldh, nnode)
!#endif
                  END DO
              END IF
           END DO
        END DO
     END IF
  END DO
  nstatold = nstatnew
  WHERE (nstatold == (-1))
     nstatold = 0
  END WHERE
  RETURN
END SUBROUTINE dbsol3
