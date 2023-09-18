 SUBROUTINE toar15(heade,lside,nrf,rhead,nelem)

 ! Create topological array
 USE ctrl_db, ONLY : npoin
 USE npo_db, ONLY : coord,label
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem !number of elements
 INTEGER (kind=4), INTENT(IN) :: nrf  !number of reinforced sides
 TYPE (ele15), POINTER :: heade   !pointer to first element in the set
 TYPE (reinf), POINTER :: rhead   !pointer to first reinforced side
 LOGICAL :: lside


 ! local variables
 LOGICAL :: found
 INTEGER (kind=4) inode,jnode,node1,node2,nodj1,nodei,nodej,iel,jel,i,j,k,ii
 INTEGER (kind=4), PARAMETER :: nextn(3) = (/2,3,1/) !cycling list
 TYPE (ele15), POINTER  :: e, ej                     !element pointers
 TYPE (reinf), POINTER :: si                         !pointer to a reinforced side
 INTEGER (kind=4) ln(6)

 !_____________________________ Find adjacent elements to reinforcement
 si => rhead            !point to first side
 ! search each side
 DO i=1,nrf             !for each side
   j = si%lnods(1)      !internal node numeration (first node of the side)
   k = si%lnods(2)      !internal node numeration (second node)
   e => heade           !point to first element
   DO iel=1,nelem           !loop
     IF(ANY(e%lnods(1:3) == j)) THEN  !if J node present
       IF(ANY(e%lnods(1:3) == k)) THEN   !if K node present
         DO ii=1,3                          !loop over each side
           IF( e%lnods(hh(3,ii)) == j .AND.  e%lnods(hh(2,ii)) == k )THEN
             si%lnods(3) = iel*10 + ii      !rigth element
           ELSE IF( e%lnods(hh(2,ii)) == j .AND.  e%lnods(hh(3,ii)) == k )THEN
             si%lnods(4) = iel*10 + ii      !left element
           ELSE
             CYCLE
           END IF
           e%si(ii)%p => si
           EXIT
         END DO
         EXIT
       END IF
     END IF
     e => e%next
   END DO
   si => si%next
 END DO
 !_____________________________

 !IF( iwrit > 0 ) &
    WRITE(lures,"(' Extended nodal connectivities for Element type RBST')",ERR=9999)
 e => heade             !point to first element
 iel = 1
 DO

   IF (.NOT.ASSOCIATED (e%next) ) EXIT !if all element processed EXIT loop

   DO inode=1,3                        !for each side in the element
     IF( e%lnods(inode+3) /= 0 )CYCLE   !if opposite already find, CYCLE
     jnode = nextn(inode)               !next local node
     node1 = e%lnods(jnode)             !first node of the side
     jnode = nextn(jnode)               !next local node
     node2 = e%lnods(jnode)             !second node of the side
     found = .FALSE.                    !initializes search
     jel = iel + 1
     ej => e%next                       !point to next element
     DO   !search an element with the same side (node1 -- node2)
       DO nodei=1,3                         !for each side
         IF( ej%lnods(nodei+3) /= 0 )CYCLE  !if side already paired, CYCLE
         nodej = nextn(nodei)            !next local node
         nodj1 = ej%lnods(nodej)         !first node of tested side
         IF(nodj1 /= node2)CYCLE         !first test, if failed, CYCLE
         nodej = nextn(nodej)            !next local node
         nodj1 = ej%lnods(nodej)         !second node of the side
         IF(nodj1 == node1) THEN         !second test
           found = .TRUE.                !side found
           e%lnods(inode+3) = ej%lnods(nodei) !opposite global node
           ej%lnods(nodei+3) = e%lnods(inode) !opposite global node
           e%lside(inode)  = jel              !adjacent element
           ej%lside(nodei) = iel              !element
           IF( ASSOCIATED(e%si(inode)%p)) THEN
             ej%si(nodei)%p => e%si(inode)%p
             IF(ej%si(nodei)%p%lnods(4) == 0 )THEN
               ej%si(nodei)%p%lnods(4) = jel*10 + nodei
             ELSE
               ej%si(nodei)%p%lnods(3) = jel*10 + nodei
             END IF
           END IF
         END IF
         EXIT                            !found or not exit loop
       END DO
       ! if twin side found or last element checked, EXIT loop
       IF (found .or. .NOT.ASSOCIATED(ej%next))EXIT
       jel = jel + 1
       ej => ej%next         !point to next element to search
     END DO
   END DO
   !IF( iwrit > 0 )THEN
      ln(1:3) = label(e%lnods(1:3))
      ln(4:6) = 0
      DO inode=1,3
        IF( ASSOCIATED(e%si(inode)%p ) )ln(inode+3) = e%si(inode)%p%numel
      END DO
      WRITE(lures,"(I7,I4,6I6)",ERR=9999) e%numel,e%matno,ln
   !END IF
   iel = iel + 1
   e => e%next               !point to next element to process
 END DO

 lside = .TRUE.              !topological array created
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE toar15
