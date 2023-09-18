 SUBROUTINE reasdf(nescv,npsdf,nesdf,ftsdf)

 !     APPLY restrictions on (translational) degrees of freedom
 USE npo_db, ONLY : ifpre,label
 USE ctrl_db, ONLY : npoin
 USE c_input
 USE nes_db
 USE outp_db, ONLY: iwrit
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN OUT) :: nescv
 INTEGER (kind=4),INTENT(OUT) :: npsdf(*),nesdf(*)
 REAL (kind=8),INTENT(OUT) :: ftsdf(*)

 INTEGER (kind=4) :: i,j,k,n,node,idof,nodo
 REAL    (kind=8) :: fact
 !      INTEGER (kind=4), PARAMETER :: nn = 1000000

 INTEGER (kind=4) chnode

   IF (iwrit == 1) WRITE(lures,"(//,5x,'Slave D.O.F. Data',//&
          &   '   Slave  DOF  Master  DOF  Factor')",ERR=9999)
   i = 0              !initializes counter for the (DATA_BASE) list
   k = 0              !initializes counter in definite list (NESDF, FTSDF)
   n = 0              !initializes number of slave DOFs  (NPSDF)
   DO
     i=i+1            !increase counter in the list
     nodo = nesdat(i)%node      !slave node label
     idof = nesdat(i)%ndof      !slave dof
     IF (iwrit == 1 ) WRITE(lures, '(i8,i5)',ERR=9999) nodo, idof
     node = chnode(nodo)        !internal node
     IF (node == 0) THEN        !??? WHAT'S THIS (FF)
       IF ( i == nnes ) EXIT    !at the end of the list, EXIT
       CYCLE                    !in the middle of the list CYCLE
     END IF
     n = n+1                    !increase number of the slave DOFs
     IF(ifpre(idof,node) < 0)THEN   !check condition
       WRITE(lures,"(' DOF already declared slave: NODE',i6,&
            &        ' DOF',i2)",ERR=9999) label(node),idof
       CALL runend('REASDF: Non consistent declaration ')
     END IF
     ifpre(idof,node) = -n      !store Slave order
     j = 0                      !initializes number of dependencies
     DO
       IF ( i == nnes ) EXIT    !if at the end of the list, EXIT
       IF ( nesdat(i+1)%factor == 0d0 ) EXIT    !if a new slave DOF exit
       i = i+1                  !increase list (DATA_BASE) counter
       k = k+1                  !increase counter in NESDF & FTSDF
       j = j+1                  !increase number of dependencies
       nodo = nesdat(i)%node    !master node
       idof = nesdat(i)%ndof    !master DOF
       fact = nesdat(i)%factor  !associated factor
       IF (iwrit == 1) WRITE(lures, '(13x,i8,i5,e12.4)',ERR=9999) &
           nodo, idof, fact          !ECHO
       node = chnode(nodo)      !internal node
       nesdf(k) = idof + 10*node     !keep value until later
       IF(ifpre(idof,node) < 0 )THEN
         WRITE(lures,"(' MASTER DOF not valid Node',i6,' dof',i2)",ERR=9999) &
                   nodo,idof
         CALL runend('REASDF: Non consistent declaration ')
       END IF
       ftsdf(k) = fact          !associated factor
     END DO
     npsdf(n+1) = npsdf(n)+j    !pointer list
     IF ( i == nnes ) EXIT      !if last value reached, EXIT (normal)

   END DO

   nescv = n                    !number of slave DOFs

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE reasdf
