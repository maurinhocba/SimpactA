 SUBROUTINE toar19(heade,lside,swapc,eside)

 ! Create topological array

 IMPLICIT NONE
 INTEGER(kind=4)       :: eside(:,:) !neighbor elements array
 TYPE (ele19), POINTER :: heade   !pointer to first element in the set
 LOGICAL :: lside,swapc

 ! local variables
 LOGICAL :: found
 INTEGER (kind=4) :: inode,jnode,node1,node2,nodj1,nodei,nodej,ns
 INTEGER (kind=4) :: elma,elmb,ielm,jelm,ie
 INTEGER (kind=4) :: i,j,k,l,m,n
 INTEGER (kind=4), PARAMETER :: nextn(3) = (/2,3,1/) !cycling list
 TYPE (ele19), POINTER  :: e, ej, ek                 !element pointers

 !_____________________________
 ielm = 1
 e => heade             !point to first element
 DO

   IF (.NOT.ASSOCIATED (e%next) ) EXIT !if all element processed EXIT loop
   DO inode=1,3                        !for each side in the element
     IF( e%lnods(inode+3) /= 0 )  CYCLE   !if opposite already find, CYCLE
     jnode = nextn(inode)               !next local node
     node1 = e%lnods(jnode)             !first node of the side
     jnode = nextn(jnode)               !next local node
     node2 = e%lnods(jnode)             !second node of the side
     found = .FALSE.                    !initializes search
     jelm = ielm + 1
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
           eside(inode,ielm) = jelm !adjacent to e  element (position)
           eside(nodei,jelm) = ielm !adjacent to ej element (position)
         END IF
         EXIT                            !found or not exit loop
       END DO
       ! if twin side found or last element checked, EXIT loop
       IF (found .OR. .NOT.ASSOCIATED(ej%next))EXIT
       jelm = jelm + 1
       ej => ej%next         !point to next element to search
     END DO
   END DO
   !WRITE(111,'(5i7)')ielm,e%numel,eside(1:3,ielm)
   ielm = ielm + 1
   e => e%next               !point to next element to process
 END DO

 ! check elements with only one neighbour
 IF( swapc )THEN
   ielm = 1
   e => heade
   DO
     IF (.NOT.ASSOCIATED (e) ) EXIT !if all element processed EXIT loop
     ns = 0
     DO inode=1,3                        !for each side in the element
       IF( e%lnods(inode+3) /= 0 )   ns = ns+1
     END DO
     IF( ns == 1) THEN         !element has one neighbour only
       inode = 1                           !initializes loop counter
       DO                                  !for each side in the element
         IF( e%lnods(inode+3) /= 0 )EXIT   !side found, EXIT
         inode = inode + 1                 !increase counter
       END DO
       i = e%lnods(inode)              !only node with opposite element
       jelm = 1
       ej => heade                    !point to first element
       loop1 : DO   !search an element with node i as the opposite node
         jnode = 1                    !initializes loop counter
         DO                            !for each side in the element
           IF( jnode > 3 )EXIT          !number of nodes exceeded, EXIT inner loop
           IF( ej%lnods(jnode+3) == i  )EXIT loop1  !node found, EXIT outer loop
           jnode = jnode + 1              !increase side counter
         END DO
         jelm = jelm + 1
         ej => ej%next                    !go to next element
       END DO loop1
       IF( ALL(ej%lnods(4:6) > 0 ))THEN   !check if element J has 3 neighbours
         node1 = nextn(inode)               !next local node in element I
         j     = e%lnods(node1)             !first node of the side I
         node1 = nextn(node1)               !next local node
         k     = e%lnods(node1)             !second node of the side I
         l     = e%lnods(inode+3)           !extra node (opposite to i)
         node2 = nextn(jnode)               !next local node on element J
         m     = ej%lnods(node2+3)          !node opposite to node k
         node2 = nextn(node2)               !next local node on element J
         n     = ej%lnods(node2+3)          !node opposite to node J
         e%lnods(1:6)  = (/ i,j,l,m,k,0 /)
         ej%lnods(1:6) = (/ k,i,l,j,n,0 /)
         elma  = eside(nextn(jnode),jelm) ! next element in patch
         elmb  = eside(node2       ,jelm) ! previous element in patch
         eside(:,ielm) = (/ elma, jelm, 0 /)
         eside(:,jelm) = (/ ielm, elmb, 0 /)
         ie = 1
         ek => heade
         loop2: DO
           IF(.NOT. ASSOCIATED(ek))EXIT
           IF    ( ie == elma )THEN
             DO inode = 4,6
               IF(ek%lnods(inode) == k) ek%lnods(inode) = i
             END DO
             DO inode = 1,3
               IF(eside(inode,ie) == jelm) eside(inode,ie) = ielm
             END DO
             elma = 0
           ELSEIF( ie == elmb )THEN
             DO inode = 4,6
               IF(ek%lnods(inode) == j) ek%lnods(inode) = i
             END DO
             elmb = 0
           END IF
           IF( (elma+elmb) == 0 ) EXIT
           ie = ie + 1
           ek => ek%next
         END DO loop2
       END IF
     END IF
     !WRITE(111,'(5i7)')ielm,e%numel,eside(1:3,ielm)
     ielm = ielm + 1
     e => e%next
   END DO
 END IF

 lside = .TRUE.              !topological array created
 RETURN
 END SUBROUTINE toar19
