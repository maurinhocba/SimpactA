 SUBROUTINE impfix (ndim1,neul1)
 !
 !   import element data set fixities
 !
 USE ctrl_db, ONLY: ndofn,neulr
 USE gvar_db, ONLY : fimpo,lab1,renum,seque,overw
 USE meshmo_db, ONLY : nodset,numpo
 USE ifx_db

 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) ::  ndim1,neul1

 INTEGER (kind=4) :: i,l,n,ns,nd,bb(7),pos(7)
 TYPE (ifx_nod), POINTER :: ifx

 READ(fimpo)ns,nd
 ! find the positions in the array
 pos(1:nd) = (/1:nd/)                !All DOFS
 IF( nd > ndim1 )THEN                !rotations are present
   IF( neulr > neul1 ) pos(nd) = nd1 !NDOFN+1
 END IF
 pos(1:nd) = pos(1:nd) + 1  !due to the label
 !              read B.C.
 DO i=1,ns                  !for each node
   READ (fimpo)n,bb(1:nd)   !read label and constraints
   IF( overw )CYCLE         !skip data
   ALLOCATE (ifx)
   ifx%ifix = 0             !initializes all values
   ! modify node label
   IF (renum)THEN
     IF( seque )THEN
       DO l=1,numpo
         IF( n == nodset(l) )THEN
           n = lab1 + l
           EXIT
         END IF
       END DO
     ELSE
       n = n + lab1
     END IF
   END IF
   ifx%ifix(1) = n                    !labels
   ifx%ifix(pos(1:nd)) = bb(1:nd)     !b.c.
   CALL add_ifx( ifx, ihead, itail )    !add node to the list
   nifx = nifx + 1
 END DO

 RETURN
 END SUBROUTINE impfix
