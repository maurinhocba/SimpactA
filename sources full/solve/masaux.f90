 SUBROUTINE masaux(nnode,lnods,nvarl,ndofe,lmass,ymass )

! auxiliar  Routine to assemble local mass into global mass

 USE npo_db, ONLY : ifpre
 USE kinc_db, ONLY : nn
 INTEGER (kind=4), INTENT(IN) :: nnode,        & !number of element nodes
                                 lnods(nnode), & !element connectivities
                                 nvarl,        & !number of element DOFs
                                 ndofe           !number of node DOFs
 REAL(kind=8), INTENT(IN) :: lmass(1)            !nvarl*(nvarl+1)/2 element matrix
 REAL(kind=8), INTENT(IN OUT) :: ymass(1)        !global matrix

 LOGICAL :: full
 INTEGER (kind=4) :: i,j,k,l,n
 INTEGER (kind=4) :: lm(nvarl)

 !to consider different number of DOFs per node
 full = nvarl == nnode*ndofe  !full matrix (same dofs at all nodes)
 l = 0                        !initializes DOF counter
 DO n = 1,nnode               !for each node
   j = lnods(n)                 !node number
   IF( j > 0)THEN               !if node exist
     DO i = 1,ndofe               !for each possible DOF
       k = ifpre(i,j)               !DOF number
       IF(full .OR. k /= 0)l = l+1    !update
       IF(k >= -nn) THEN
         IF(full .OR. k /= 0) lm(l) = k
       ELSE
         lm(l) = 0
       END IF
     END DO
   ELSE IF(full)THEN
     lm(l+1:l+ndofe) = 0
     l = l+ndofe
   ELSE  !ONLY for TTTL
     lm(l+1:l+3) = 0
     l = l+3
   END IF
 END DO
 CALL ensmat(nvarl,lm(1),lmass(1),ymass(1))

 RETURN
 END SUBROUTINE masaux
