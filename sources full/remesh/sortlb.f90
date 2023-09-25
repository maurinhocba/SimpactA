       SUBROUTINE sortlb(nn,iarr,indx)
       !***************************************************************
       ! This subroutine sorting the array and
       !***************************************************************
       IMPLICIT NONE
       INTEGER(kind=4),INTENT(IN)  :: nn
       INTEGER(kind=4),INTENT(IN)  :: iarr(nn)
       INTEGER(kind=4),INTENT(OUT) :: indx(nn)
       INTEGER(kind=4) :: a
       INTEGER(kind=4) :: i,j,k,indext,jstack,l,r
       INTEGER(kind=4) :: istack(nn)
       DO i=1,nn    ! create the indx array WC
         indx(i) = i
       END DO
       jstack = 0
       l = 1
       r = nn
       DO
         IF(r-l < NN)THEN
           DO j=l+1,r
             indext = indx(j)
             a      = iarr(indext)
             DO i=j-1,1,-1
               IF(iarr(indx(i)) <= a)EXIT
               indx(i+1) = indx(i)
             END DO
             indx(i+1) = indext
           END DO
           IF(jstack == 0)RETURN
           r = istack(jstack)
           l = istack(jstack-1)
           jstack = jstack-2
         ELSE
           k = (l+r)/2
           CALL swap(indx(k),indx(l+1))
           CALL icomp_xchg(indx(l),indx(r))
           CALL icomp_xchg(indx(l+1),indx(r))
           CALL icomp_xchg(indx(l),indx(l+1))
           i = l+1
           j = r
           indext = indx(l+1)
           a = iarr(indext)
           DO
             DO
               i = i+1
               IF(iarr(indx(i)) >= a)EXIT
             END DO
             DO
               j = j-1
               IF(iarr(indx(j)) <= a)EXIT
             END DO
             IF(j < i)EXIT
             CALL swap(indx(i),indx(j))
           END DO
           indx(l+1) = indx(j)
           indx(j)   = indext
           jstack = jstack+2
           IF(jstack > nn) CALL runen3('INDEX.F90: stack problem in sorting')
           IF(r-i+1>=j-l)THEN
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
           ELSE
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l = i
           END IF
         END IF
       END DO

       CONTAINS

         SUBROUTINE icomp_xchg(i,j)
         ! Swap component in indexing array
         INTEGER(kind=4),INTENT(INOUT) :: i,j
         INTEGER(kind=4)               :: swp
         IF(iarr(j) < iarr(i))THEN
           swp = i
           i   = j
           j   = swp
         END IF
         END SUBROUTINE icomp_xchg

         SUBROUTINE swap(a,b)
         ! Swap the contents of a and b.
         INTEGER(kind=4),INTENT(INOUT) :: a,b
         INTEGER(kind=4)               :: dum
         dum = a
         a   = b
         b   = dum
         END SUBROUTINE swap

       END SUBROUTINE sortlb
