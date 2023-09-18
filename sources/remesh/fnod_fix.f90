      SUBROUTINE fnod_fix(nnodbd,nodcre)
!--------------------------------------------------------------------
!     get ifx conditions for new boundary nodes in a refined set
!--------------------------------------------------------------------     
      USE ctrl_db, ONLY: ndofn, ndime
      USE outp_db, ONLY: iwrit    
      USE c_input  
      USE ifx_db,  ONLY: nd1, nifx, ifx_nod, ihead, itail, srch_ifx, add_ifx
      USE npo_db,  ONLY: oldlb      
      USE meshmo_db,ONLY: maxnn      
      IMPLICIT NONE
      !--- Dummy variables
      INTEGER(kind=4),POINTER:: nodcre(:,:), & !neighbors nodes between which a new node is generated
                                nnodbd(:)      !nnodbd = 0 If created node is not a boundary node
                                               !       = 1 If created node is a boundary node
      !--- Local variables
      INTEGER(kind=4):: ijfix(nd1,2)   ! fixity in i-j segment (2 nodes segments)                        
      INTEGER(kind=4):: i, j, n, dim1, inode, nfix, newno
      INTEGER(kind=4),POINTER:: ifx_b(:,:) !(nd1,"new boundary nodes") boundary conditions
      LOGICAL:: found, newfx
      TYPE (ifx_nod), POINTER :: ifx, anter, posic

      INTERFACE
        INCLUDE 'inc_mint.h'
      END INTERFACE

      nfix  = 0 ! init new fixities conditions
      DO i=1,SIZE(nnodbd)  ! search for new boundary nodes          
         IF( nnodbd(i) == 0 )CYCLE ! if node isn't in boundary CYCLE
         dim1 = 0 ! initialize
         newno = maxnn + i ! label of the new node
         !extract i-j neighbors nodes (segment) and check boundary conditions
         ijfix = 0 !initialize
         DO inode=1,2 !search in fix_db for "i" an "j" segment labels nodes
           n = oldlb(nodcre(inode,i)) !search for old label of segment nodes
           CALL srch_ifx(ihead, anter, posic, n, found)
           IF (found) ijfix(1:nd1,inode) = posic%ifix(2:nd1+1) !assign nodal boundary cond.
         END DO

         newfx = .FALSE. !initialize
         DO j=1,nd1 ! check if ALL segment is FIXED in some dof
            IF( ALL(ijfix(j,1:2) /= 0) )newfx = .TRUE.
         END DO

         IF(newfx)THEN
           nfix = nfix + 1 ! increase number of new boundary nodes conditions
           IF( nfix == 1 )dim1 = nd1+1 !for allocation ifx_b only
           CALL inc_mint(ifx_b,dim1,1,dflt=0) !get memory for each new fixity conditions
           ifx_b(1,nfix) = newno ! label of new node with boundary condition
           DO j=1,nd1 ! only if i-j segment is FIXED then save this condition
             IF( ALL(ijfix(j,1:2) /= 0) )ifx_b(2:nd1+1,nfix) = ijfix(:,1)
           END DO
         END IF

      END DO

      ! add new boundary nodes to fixation database
      IF (iwrit == 1) THEN
        WRITE(lures,"(/,' Added boundary Conditions',/)",ERR=9999)      
        IF(ndime == 2) WRITE(lures,"(/,'     Node   XYA')",ERR=9999)
        IF(ndime == 3) WRITE(lures,"(/,'     Node   XYZABG')",ERR=9999)
      END IF
      
      nifx = nifx + nfix !increment total number of nodes in the fixity list
      DO i=1,nfix
        ALLOCATE (ifx)                !get memory      
        ifx%ifix(1) = ifx_b(1,i)          !node label
        ifx%ifix(2:nd1+1)= ifx_b(2:nd1+1,i)  !restrictions
        IF (iwrit == 1) WRITE(lures,'(1i9,3x,7i1)',ERR=9999) ifx%ifix(1:nd1+1)  !echo
        CALL add_ifx( ifx, ihead, itail )  !add to list
      END DO

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE fnod_fix



