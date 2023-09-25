 SUBROUTINE rigbdy ( nmast,slaves,nnode)
 !
 !     This routine computes:
 !          center of mass and mass of a rigid body
 !          inertia tensor, principal (eigen-) values a directions
 !          nodal dependencies
 !
 USE ctrl_db, ONLY : ndime,neulr
 USE npo_db,  ONLY : emass,coora,coord,label,euler,naeul
 USE outp_db, ONLY: iwrit
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: &
                   nnode,             &! number of nodes
                   slaves(:),         &! list of nodes
                   nmast               ! master node

 INTEGER (kind=4) n1,j,n,k
 REAL (kind=8) xc(ndime),d(ndime),mm
 REAL (kind=8),ALLOCATABLE :: mas(:),Iner(:,:)
 LOGICAL first,eul

 IF(nnode < 1)RETURN         !Exit if no elements
 k = 1                       !inertia tensor dimension
 IF(ndime == 3) k = ndime    ! only Z for 2-D, X-Y-Z for 3-D
 eul = neulr > 0
 IF( eul ) eul = naeul(nmast)
 !      Mass value and position  Inertia tensor,  nmast=number of master node
 ALLOCATE( mas(ndime+1) , iner(k,k) )
 mas  = 0.0d0     !initializes
 iner = 0.0d0

 !     First loop on elements to compute center of mass
 DO n = 1,nnode                     !For each node
   n1 = slaves(n)                   !node
   DO j=1,ndime                     !for each direction
     mas(j) = mas(j) + coora(j,n1)*emass(j,n1)   !mass 1st moment
   END DO
   mas(1+ndime) = mas(1+ndime) + emass(1,n1) !mass
 END DO

 !     loop on nodes for the computation of center of mass

 d = coora(:,nmast) - coord(:,nmast)
 first = SQRT(DOT_PRODUCT(d,d)) < 1e-12   !displ are zero
 IF(mas(1+ndime) > 0)THEN  !If the body has mass
   DO j=1,ndime                          !for each direction
     coora(j,nmast) = mas(j)/mas(1+ndime)   !center of mass
     IF( first ) coord(j,nmast) = coora(j,nmast)
   END DO
   IF(iwrit == 1) &
    WRITE(lures,"(//4x,'Rigid_B_Node',5x,'X',9x,'Y',9x,'Z',6x,'Mass',/, &
                & i16,4f10.3)",ERR=9999) label(nmast),coora(:,nmast),mas(ndime+1)

 ELSE
   DEALLOCATE( mas , iner )
   IF(emass(1,nmast) > 0d0) RETURN      !exit, mass is concentrated at master node
   WRITE(lures,"('Rigid body without mass, Master node:',i10)",ERR=9999) label(nmast)
   CALL runen3('RIGBDY: a rigid body without mass defined')
 END IF
 !     Second loop on elements to compute rotational inertia
 IF( eul )THEN
   xc = coora(1:ndime,nmast)         !nodal coordinates of center of mass
   DO n = 1,nnode              !for each node in the element
     n1 = slaves(n)            !node number
     d = coora(1:ndime,n1) - xc    !distance to center of mass
     mm = emass(1,n1)          !nodal mass
     IF(ndime == 2)THEN        !2-d problems, only Z inertia
       iner(1,1) = iner(1,1) + mm*(d(1)*d(1)+d(2)*d(2)) + emass(3,n1)
     ELSE                      !3-d problems, whole inertia tensor (3x3)
       iner(1,1) = iner(1,1) + mm*(d(2)*d(2)+d(3)*d(3)) + emass(4,n1)
       iner(1,2) = iner(1,2) - mm*d(1)*d(2)
       iner(1,3) = iner(1,3) - mm*d(1)*d(3)
       iner(2,1) = iner(2,1) - mm*d(2)*d(1)
       iner(2,2) = iner(2,2) + mm*(d(1)*d(1)+d(3)*d(3)) + emass(5,n1)
       iner(2,3) = iner(2,3) - mm*d(2)*d(3)
       iner(3,1) = iner(3,1) - mm*d(3)*d(1)
       iner(3,2) = iner(3,2) - mm*d(3)*d(2)
       iner(3,3) = iner(3,3) + mm*(d(1)*d(1)+d(2)*d(2)) + emass(6,n1)
     END IF
   END DO

   !     loop on r-b for the computation of principal direction of inertia tensor

   IF(iwrit == 1)WRITE(lures,"('Rigid_b_node inertial tensors')",ERR=9999)
   IF(ndime == 2)THEN   !for 2-d problems, Z is the principal axis
     IF(iwrit == 1)WRITE(lures,"(i10,e15.5)",ERR=9999) label(nmast),iner(1,1)
   ELSE                 !for 3-d problems, compute directions
     CALL princi(ndime,iner(1,1),euler(1,nmast))  ! ==> euler(:,nmast)
     IF(iwrit == 1) WRITE(lures,"(i10,3e15.5,/,(10x,3f15.5))",ERR=9999) label(nmast), &
           &    (iner(j,j),j=1,3),euler(1:neulr,nmast)
   END IF
 END IF
 DEALLOCATE( mas , iner )

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rigbdy
