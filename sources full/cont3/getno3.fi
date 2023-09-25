 SUBROUTINE getno3(nsegm,nnseg,lcseg,ncnoo,lcnoo,ncnod,lcnod,newn,ends)

 !.... set list of nodes from connectivities for a surface
 !....     compare with existing list
 !
 !.... input
 !       lcnoo(ncnoo) : nodes in the original slave surface
 !       ncnoo = number of nodes in the original surface
 !....   nsegm = number of segments of the new  surface considered
 !....   nnseg = number of nodes per segment of the surface considered
 !....   lcseg(inseg,icseg) = global node number for the local element node
 !....                        [inseg] of the segment [icseg]
 !.... output
 !....   ncnod = number of nodes of the surface considered
 !....   lcnod(ncnod) = global node number for the local node [icnod] of the
 !....                  surface considered
 !....   newn = number of new nodes
 !....   ends(optional) = conections between new nodes and old nodes

 IMPLICIT NONE
 !     arguments
 INTEGER(kind=4),INTENT(IN) :: nsegm,nnseg,lcseg(:,:),ncnoo
 INTEGER(kind=4),INTENT(IN OUT) :: lcnoo(:)
 INTEGER(kind=4),INTENT(OUT) :: ncnod,lcnod(:),newn
 INTEGER (KIND=4), OPTIONAL  :: ends(:,:)
 !     local variables
 INTEGER (kind=4)  jm,icseg,inseg,i,j,k,im
 LOGICAL, ALLOCATABLE :: exi(:)
 LOGICAL :: finde

 ALLOCATE ( exi(ncnoo) )
 exi  = .FALSE.
 finde = PRESENT(ends)
 IF( finde ) ends = 0

 DO icseg = 1,nsegm            !for each segment in the mesh
   dos : DO inseg = 1,nnseg    !for each node in the segment
     jm = lcseg(inseg,icseg)   !node
     DO i=1,ncnoo              !search node in the old list
       IF(jm == lcnoo(i))THEN  !if node exists in the old list
         exi(i) = .TRUE.       !set as existing
         CYCLE dos             !exit loop
       END IF
     END DO
     !search in the new list
     i = 0                     !initializes counter
     DO
       i = i+1                 !increments counter
       IF( i > ncnod ) THEN    !check
         ncnod = i             !if not counted, increment counter
         lcnod(i) = jm         !add to list
         EXIT                  !exit loop
       END IF
       IF( lcnod(i) == jm ) EXIT  !if node already in the list, exit loop
     END DO
     ! keep conexions between new nodes and old nodes
     IF( finde ) THEN                !if array ENDS present
       DO j=1,nnseg                  !for each node in the segment
         IF( j == inseg )CYCLE       !not for the new node of course
         im = lcseg(j,icseg)         !a vertex node
         IF( ANY( lcnoo == im ) )THEN  !if node exists in the old mexh
           DO k=1,4                    !for each possibilitie
             IF(ends(k,i) == 0 )THEN   !if not assigned
               ends(k,i) = im          !assign
               EXIT                    !exit loop
             ELSE
               IF(ends(k,i) == im)EXIT !if node already kept
             END IF
           END DO
         END IF
       END DO
     END IF
   END DO dos
 END DO
 ! add old nodes to new nodes
 newn = ncnod !keep number of new nodes
 DO j=1,ncnoo !for each old node
   IF( exi(j) )THEN             !if kept in new mesh
     ncnod = ncnod + 1          !increment counter
     lcnod(ncnod) = lcnoo(j)    !keep node number
     lcnoo(j) = ncnod           !remember position in new list
   ELSE                         !for a deleted node
     lcnoo(j) = 0               !set 0
   END IF
 END DO

 RETURN
 END SUBROUTINE getno3

