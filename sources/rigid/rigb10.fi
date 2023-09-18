 SUBROUTINE rigb10 (nelem,lnods,nmast,nnode)
 !
 !     This routine computes:
 !          center of mass and mass of each rigid body
 !          inertia tensor, principal (eigen-) values & directions
 !
 USE ctrl_db, ONLY : ndime,neulr
 USE npo_db,  ONLY : naeul,euler,coora,coord,emass
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) ::   &
      nelem,        &     ! Number of elements
      nnode,        &     ! number of nodes per element
      lnods(:,:),   &     ! Connectivities
      nmast               ! Master node

 INTEGER (kind=4) iel,n1,j,n,k,im,jm
 REAL (kind=8) xc(ndime),d(ndime),mm
 REAL (kind=8),ALLOCATABLE :: mas(:),Iner(:,:)
 LOGICAL :: first,eul
 LOGICAL, ALLOCATABLE :: used(:)

 IF(nelem < 1)RETURN         !Exit if no elements
 k = 1                       !inertia tensor dimension, only Z for 2 D
 IF(ndime == 3) k = ndime    !X-Y-Z for 3-D
 !      Mass value and position  Inertia tensor,  nrbod=number of R-B
 im = MINVAL(lnods)           !minimum node number
 jm = MAXVAL(lnods)           !maximum node number
 ALLOCATE( mas(ndime+1) , iner(k,k), used(im:jm) )
 mas  = 0.0d0     !initializes
 iner = 0.0d0
 used = .FALSE.

 !     First loop on elements to compute center of mass
 DO iel = 1,nelem                   !For each element
   DO n = 1,nnode                   !for each node in the element
     n1 = lnods(n,iel)              !node
     IF( used(n1) )CYCLE            !if already considered
     DO j=1,ndime                   !for each direction
       mas(j) = mas(j) + coora(j,n1)*emass(j,n1)   !mass 1st moment
     END DO
     mas(1+ndime) = mas(1+ndime) + emass(1,n1) !mass
     used(n1) = .TRUE.
   END DO
 END DO

 !     loop on elements for the computation of center of mass
 WRITE(lures,"(//4x,'Node',7x,'X',14x,'Y',14x,'Z',12x,'Mass')",ERR=9999)

 ! ADD mass in the Rigid Bodie node
 DO j=1,ndime                   !for each direction
   mas(j) = mas(j) + coora(j,nmast)*emass(j,nmast)   !mass 1st moment
 END DO
 mas(ndime+1) = mas(ndime+1) + emass(1,nmast)
 d = coora(:,nmast) - coord(:,nmast) !displacement
 first = SQRT(DOT_PRODUCT(d,d)) < 1e-12
 IF(mas(1+ndime) > 0)THEN  !If the body has mass
   DO j=1,ndime                !for each direction
     coora(j,nmast) = mas(j)/mas(1+ndime) !center of mass
   END DO
   IF( first)coord(:,nmast) = coora(:,nmast)
   WRITE(lures,"(i8,4f15.10)",ERR=9999) nmast,coora(:,nmast),mas(ndime+1)
 END IF

 eul = neulr > 0
 IF( eul ) eul = naeul(nmast)
 IF( eul ) THEN !IF local system is defined for the master node
   used = .FALSE.       !Initializes
   !     Second loop on elements to compute rotational inertia
   DO iel = 1,nelem              !for each element
     xc = coora(1:ndime,nmast)         !nodal coordinates of center of mass
     DO n = 1,nnode              !for each node in the element
       n1 = lnods(n,iel)         !node number
       IF( used(n1) )CYCLE       !if already considered, next node
       d = coora(1:ndime,n1) - xc    !distance to center of mass
       mm = emass(1,n1)          !nodal mass
       IF(ndime == 2)THEN        !2-d problems, only Z inertia
         iner(1,1) = iner(1,1) + mm*(d(1)*d(1)+d(2)*d(2))
       ELSE                      !3-d problems, whole inertia tensor (3coora3)
         iner(1,1) = iner(1,1) + mm*(d(2)*d(2)+d(3)*d(3))
         iner(1,2) = iner(1,2) - mm*d(1)*d(2)
         iner(1,3) = iner(1,3) - mm*d(1)*d(3)
         iner(2,1) = iner(2,1) - mm*d(2)*d(1)
         iner(2,2) = iner(2,2) + mm*(d(1)*d(1)+d(3)*d(3))
         iner(2,3) = iner(2,3) - mm*d(2)*d(3)
         iner(3,1) = iner(3,1) - mm*d(3)*d(1)
         iner(3,2) = iner(3,2) - mm*d(3)*d(2)
         iner(3,3) = iner(3,3) + mm*(d(1)*d(1)+d(2)*d(2))
       END IF
       used(n1) = .TRUE.         !global node is part of the rigid body
     END DO
   END DO

   !     loop on r-b for the computation of principal direction of inertia tensor

   WRITE(lures,"('    node inertial tensors')",ERR=9999)
   d = coord(1:ndime,nmast) - coora(1:ndime,nmast)    !distance to center of mass
   mm = emass(1,nmast)          !nodal mass
   IF(ndime == 2)THEN   !for 2-d problems, Z is the principal axis
     iner(1,1) = iner(1,1) + mm*(d(1)*d(1)+d(2)*d(2))
     WRITE(lures,"(i7,e15.5)",ERR=9999) nmast,iner(1,1)
   ELSE                 !for 3-d problems, compute directions
     iner(1,1) = iner(1,1) + mm*(d(2)*d(2)+d(3)*d(3))
     iner(1,2) = iner(1,2) - mm*d(1)*d(2)
     iner(1,3) = iner(1,3) - mm*d(1)*d(3)
     iner(2,1) = iner(2,1) - mm*d(2)*d(1)
     iner(2,2) = iner(2,2) + mm*(d(1)*d(1)+d(3)*d(3))
     iner(2,3) = iner(2,3) - mm*d(2)*d(3)
     iner(3,1) = iner(3,1) - mm*d(3)*d(1)
     iner(3,2) = iner(3,2) - mm*d(3)*d(2)
     iner(3,3) = iner(3,3) + mm*(d(1)*d(1)+d(2)*d(2))
     CALL princi(ndime,iner(1,1),euler(1,nmast))  ! ==> euler(:,nmast)
     WRITE(lures,"(i7,3e15.5,/,(10x,3f12.7))",ERR=9999) nmast,            &
          &     (iner(j,j),j=1,3),euler(1:neulr,nmast)
   END IF

 END IF
 DEALLOCATE( mas , iner, used )  !release auxiliar memory
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rigb10
