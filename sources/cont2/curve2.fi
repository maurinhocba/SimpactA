 SUBROUTINE curve2(x,surf,np)
 !
 ! Compute nodal normals and segment Z-increments of a surface
 !
 IMPLICIT NONE
 INTEGER (KIND=4), INTENT (IN) :: np
 REAL (KIND=8), POINTER :: x(:,:)
 TYPE (surf2_db) :: surf
 ! INTENT(OUT) :: surf%cu(nsegm) :segment curvatures for Master-Curved surfaces
 !                surf%tn(2,ncnod) :nodal normals for Slave-Curved surfaces or self contact
 !                surf%tc(2,nsegm) :segment normals for self contact
 !                surf%xc(2,nsegm) :segment position of center
 !local variables
 INTEGER (KIND=4), PARAMETER :: nxn(2) = (/ 2,1 /)
 INTEGER (KIND=4) :: ncnod,nsegm,i,j,k,n,nn(2)
 INTEGER (KIND=4), POINTER :: nr(:)
 REAL (KIND=8) :: a(2),t(2),long,tt(2,2)

 LOGICAL :: master,slave,auto,curved
 REAL (KIND=8), POINTER :: tn(:,:),tc(:,:)


 ncnod = surf%ncnod            !number of nodes in the surface
 nsegm = surf%nsegm            !number of segments defining surface

 slave  = surf%iscod      !If acts as slave
 master = surf%imcod      !If acts as master
 auto   = surf%auto       !If self contact
 curved = surf%curved     !If surface treated as curved

 IF(master) THEN
   IF(curved) THEN
     IF( .NOT.ASSOCIATED (surf%cu ) )ALLOCATE(surf%cu(nsegm))
   END IF
   !IF( .NOT.ALLOCATED (surf%xc ) )ALLOCATE(surf%xc(2,nsegm))
   IF( auto )THEN
     IF( .NOT.ASSOCIATED (surf%tc ) )ALLOCATE(surf%tc(2,nsegm))
     tc => surf%tc
   ELSE
     ALLOCATE( tc(2,nsegm) )
   END IF
 END IF
 !IF( .NOT.ASSOCIATED (surf%nr ) )ALLOCATE(surf%nr(np))
 !nr => surf%nr
 ALLOCATE( nr(np) )
 ALLOCATE( tn(2,ncnod) )
 tn = 0d0                      !average normals initialization
 nr = 0                        !inverse relation initialization

 DO i=1,ncnod                  !generate inverse relation (INEFFICIENT)
   ! it stores at the global position the local position in LCNOD
   nr(surf%lcnod(i)) = i
 END DO
 ! Compute average normals
 DO n=1,nsegm                  !for each segment
   nn(:) = surf%lcseg(:,n)     !connectivities of the segment
   t(1) = -x(2,nn(2)) + x(2,nn(1))  !normal*l
   t(2) =  x(1,nn(2)) - x(1,nn(1))  !normal*l
   long = DOT_PRODUCT(t,t)          !squared length

   DO i=1,2                         !sums on average normal weighted by the inverse of the area
     k = nr(nn(i))                  !local numeration of the node
     tn(:,k) = tn(:,k) + t/long     !sums
   END DO
   IF(master)THEN
     tc(:,n) = t              !store normals Weigthed by the length
     surf%xc(:,i) = (x(:,nn(1)) + x(:,nn(2)) )/2d0 !compute centers
   END IF
 END DO
 ! compute unitary normals
 DO n=1,ncnod
   t = tn(:,n)                       !weighted normal at nodes
   long = SQRT(DOT_PRODUCT(t,t))     !length
   tn(:,n) = t/long                  !unit vector
 END DO

 IF( master)THEN !If surface acts as MASTER compute normal increments at the side
   DO n=1,nsegm                      !for each segment
     t = tc(:,n)                     !restore normal
     long = SQRT(DOT_PRODUCT(t,t))   !compute length
     t = t/long                      !unit vector
     nn(:) = surf%lcseg(:,n)         !connectivities of the segment
     DO i=1,2                        !for each node in the segment
       k = nr(nn(i))                 !local numeration of the node
       tt(:,i) = tn(:,k)             !nodal averaged normal
       !COS of angle between nodal average normal and segment normal
       a(i) = DOT_PRODUCT(tt(:,i),t)
     END DO
     !         Compute Z increments  along sides
     IF(curved) THEN
       surf%cu(n) = DOT_PRODUCT(tt(:,2)-tt(:,1),(/ tc(2,n), -tc(1,n) /) )
       k = 0                           !initializes
       DO i=1,2
         !modify Z-incrementes for plane segments
         IF(ABS(a(i) - 1d0) < 1D-5)THEN      !same normal for node and segment
           j = nxn(i)                        !side nodes
           !if a side element exist
           IF(surf%nhseg(i,n) > 0 .OR. surf%nhseg(j,n) > 0) surf%cu(n) = 0d0
           EXIT
         END IF
       END DO
     END IF
     tc(:,n) = t                     !keep unit normal
   END DO
   IF(master .AND. .NOT.auto) DEALLOCATE (tc) !release memory
 END IF

 IF( (slave .AND. curved) .OR. auto  )THEN
   surf%tn => tn     !keep normals if the surface acts as SLAVE
 ELSE
   DEALLOCATE( tn, surf%lcnod ) !release memory
 END IF

 DEALLOCATE (nr)     !release temporary memory

 RETURN
 END SUBROUTINE curve2
