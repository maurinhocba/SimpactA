 SUBROUTINE toar14(lside,nelem,heade)

 ! Create topological array

 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem !number of elements
 TYPE (ele14), POINTER :: heade   !pointer to first element in the set
 LOGICAL :: lside

 ! local variables
 INTEGER (kind=4) maxsd,i,ii,j,jj,k,ll,l,n,nsd,j1,k1,ielem,jelem,nn,jjj,nbs,nl(6)
 INTEGER (kind=4), PARAMETER :: nextn(3) = (/2,3,1/) !cycling list
 TYPE (ele14), POINTER  :: e, ej                     !element pointers
 LOGICAL :: first,err
 INTEGER (kind=4), ALLOCATABLE :: ns(:,:),ln(:,:),lock(:)

 !_____________________________
 maxsd = 2*nelem+1      !maximum number of sides
 ALLOCATE (ns(7,maxsd)) !auxiliar array (addmits 4 elements over one side)
 ns = 0                 !initializes array
 nsd = 0                !initializes number of sides
 e => heade             !point to first element
 ! count appearence of each side
 DO i=1,nelem           !for each element
   DO ii=1,3            !for each side
     j = e%lnods(kk(3,ii))     !internal node numeration (first node of the side)
     k = e%lnods(kk(2,ii))     !internal node numeration (second node)
     first = j < k             !J is the node with lower numeration
     l = MIN(j,k)              !use the lower value
     k = MAX(j,k)              !greater value
     j = l                     !lower value
     jj = 1                    !initializes pointer
     DO                    !search in existing sides
       IF( jj > nsd )THEN  !all sides checked ?
         nsd = jj          !increase number of sides
         ns(1,jj) = j      !keep side nodes
         ns(2,jj) = k
         ns(3,jj) = 4      !position
         EXIT
       ELSE IF(ns(1,jj) == j .AND. ns(2,jj) == k )THEN  !side exists
         ns(3,jj) = ns(3,jj) + 1   !position
         EXIT
       END IF
       jj = jj + 1         !next side
     END DO
     ! keep element and side in one value
     l = ns(3,jj)          !position in array NS
     IF( first )THEN   !ordered
       ns(l,jj) =   ii+10*i    !keeps ii = side in element and i = element
     ELSE              !reverse order
       ns(l,jj) = -(ii+10*i)
     END IF
   END DO
   e => e%next          !go to next element
 END DO
 ! generates data base
 e => heade             !point to first element
 ielem = 1              !initializes element I
 ej => heade            !point to first element
 jelem = 1              !initializes element J
 nbs = 0                !initializes number of branching nodes
 err = .FALSE.
 DO i=1,nsd             !for each side
   nn = ns(3,i)-3       !number of connected elements
   IF( nn > 2) THEN     !if a branching side
     nbs = nbs + 1        !increase number

     IF(nbs == 1)THEN   !for first branching side
       WRITE(lures,"(' error in connectivities, more than 2 elements sharing a side (Branching)')",ERR=9999)
       WRITE(*    ,"(' error in connectivities, more than 2 elements sharing a side (Branching)')")
       WRITE(55   ,"(' error in connectivities, more than 2 elements sharing a side (Branching)')",ERR=9999)
     END IF
     WRITE(lures,"(' Side nodes: ',2i7)",ERR=9999) ns(1:2,i)    !side nodes j-k (ordered)
     DO ii=1,nn                     !for each element on this side
       jj = ns(ii+3,i)     !element + side
       j = MOD(ABS(jj),10) !side
       k = ABS(jj-j)/10    !element
       WRITE(lures,"(i7,i3)",ERR=9999) k,j  !element and side
     END DO

   ELSE IF( nn == 2 )THEN          !standard element

     !  first element
     jj = ns(4,i)          !element + side
     j = MOD(ABS(jj),10)   !side
     k = ABS(jj-j)/10      !element
     jjj= jj/k           !+/-1
       !  point to element k
     IF( k > ielem )THEN
       n = ielem
     ELSE
       e => heade
       n = 1
     END IF
     DO ielem=n,k-1
       e => e%next
     END DO
     ! second element
     ll = ns(5,i)          !element + side
     j1 = MOD(ABS(ll),10)  !side1
     k1 = ABS(ll-j1)/10    !element1
       !  point to element k1
     IF( k1 > jelem )THEN
       n = jelem
     ELSE
       ej => heade
       n = 1
     END IF
     DO jelem=n,k1-1
       ej => ej%next
     END DO
     !  store information
     IF( jjj*ll > 0 )THEN
       IF( .NOT.err )THEN
         WRITE(lures,"('  WARNING, non-orientated surface found ')",ERR=9999)
         WRITE(*    ,"('  WARNING, non-orientated surface found ')")
         WRITE(55   ,"('  WARNING, non-orientated surface found ')",ERR=9999)
       END IF
       WRITE(lures,"('  Elements  ',i6,' and ',i6)",ERR=9999) e%numel,ej%numel
       e%lside(j)   = -jelem   !keep side elements
       ej%lside(j1) = -ielem   !keep side elements
       err = .TRUE.
     ELSE
       e%lside(j)  = jelem    !keep side elements
       ej%lside(j1) = ielem   !keep side elements
     END IF
     e%lnods(j+3)   = ej%lnods(j1) !opposite global node
     ej%lnods(j1+3) = e%lnods(j)   !opposite global node
   END IF
 END DO
 IF( nbs > 0) STOP 'branched surface not allowed'
 DEALLOCATE (ns)
 IF( err ) THEN !reverse connectivities found
   ALLOCATE(ln(3,nelem),lock(nelem))
   lock = 0
   e => heade
   DO i=1,nelem
     ln(1:3,i) = e%lside(1:3)
     e => e%next
   END DO
   CALL lock14(ln,1,lock,1)
   e => heade
   DO i=1,nelem
     SELECT CASE (lock(i))
     CASE (-1 ) !reverse things
       nl = e%lnods(1:6)
       e%lnods(1:6) = nl((/ 1,3,2,4,6,5  /))
       e%lside(1) = ABS(ln(1,i))
       e%lside(2) = ABS(ln(3,i))
       e%lside(3) = ABS(ln(2,i))
     CASE (0)
       STOP 'element not accessed'
     CASE (1 ) !check
       DO j=1,3
         IF( e%lside(j) < 0 ) e%lside(j) = -e%lside(j)
       END DO
     END SELECT
     e => e%next
   END DO
   DEALLOCATE (ln,lock)
 END IF
 lside = .TRUE.              !topological array created
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE toar14
   RECURSIVE SUBROUTINE lock14(ln,i,lock,signo)
   !
   IMPLICIT none

     INTEGER (kind=4),INTENT(IN) ::  ln(:,:),i,signo
     INTEGER (kind=4),INTENT(IN OUT) :: lock(:)
     INTEGER (kind=4) :: j,k

     DO j=1,3           !for each side
       k = ln(j,i)      !side element
       IF( k == 0 )CYCLE   !no side element cycle
       IF( k > 0 )THEN        !side element and coincident normals
         IF( lock(k) == 0 )THEN     !not locked
           lock(k) = signo            !lock
           CALL lock14(ln,k,lock,signo) !and continue
         ELSE  !alreadky locked
           CYCLE
         END IF
       ELSE
         IF( lock(-k) == 0 )THEN     !not locked
           lock(-k) = -signo
           CALL lock14(ln,-k,lock,-signo) !and continue
         ELSE  !alreadky locked
           CYCLE
         END IF
       END IF
     END DO

   RETURN
   END SUBROUTINE lock14
