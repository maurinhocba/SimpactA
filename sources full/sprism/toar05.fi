 SUBROUTINE toar05(lface,nelem,heade)

 ! Create topological array

 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem !number of elements
 TYPE (ele05), POINTER :: heade   !pointer to first element in the set
 LOGICAL :: lface

 ! local variables
 INTEGER (kind=4) maxsd,i,ii,j,jj,k,ll,l,n,nsd,j1,k1,ielem,jelem,nn,lf(4),ls(4),kk(1)
 INTEGER (kind=4), PARAMETER :: face(4,3) = (/ 2,3,6,5, 3,1,4,6, 1,2,5,4 /) !
 TYPE (ele05), POINTER  :: e, ej                     !element pointers
 INTEGER (kind=4), ALLOCATABLE :: ns(:,:)

 !_____________________________
 maxsd = 2*nelem+1      !maximum number of faces
 ALLOCATE (ns(7,maxsd)) !auxiliar array (addmits 2 elements over each face)
 ns = 0                 !initializes array
 nsd = 0                !initializes number of faces
 e => heade             !point to first element
 ! count appearence of each face
 DO i=1,nelem           !for each element
   DO ii=1,3            !for each face
     lf = e%lnods(face(:,ii))
     DO j = 4,1,-1
       kk = MAXLOC(lf)
       ls(j) = lf(kk(1))
       lf(kk(1)) = 0
     END DO
     jj = 1                !initializes pointer
     DO                    !search in existing faces
       IF( jj > nsd )THEN  !all faces checked ?
         nsd = jj          !increase number of faces
         ns(1:4,jj) = ls   !keep face nodes
         ns(5,jj) = 6      !position
         EXIT
       ELSE IF(ns(1,jj) == ls(1) .AND. ns(2,jj) == ls(2) .AND. ns(3,jj) == ls(3) .AND. ns(4,jj) == ls(4) )THEN  !face exists
         ns(5,jj) = ns(5,jj) + 1   !position
         EXIT
       END IF
       jj = jj + 1         !next face
     END DO
     ! keep element and face in one value
     l = ns(5,jj)          !position in array NS
     ns(l,jj) = ii+10*i    !keeps ii = face in element and i = element
   END DO
   e => e%next          !go to next element
 END DO
 ! generates data base
 e => heade             !point to first element
 ielem = 1              !initializes element I
 ej => heade            !point to first element
 jelem = 1              !initializes element J
 DO i=1,nsd             !for each face
   nn = ns(5,i)-5       !number of connected elements
   IF( nn == 2 )THEN          !standard element
     !  first element
     jj = ns(6,i)          !element + face
     j = MOD(jj,10)     !face
     k = (jj-j)/10      !element
     !  point to element k
     IF( k > ielem )THEN
       n = ielem
     ELSE
       e => heade
       n = 1
     END IF
     DO ielem=n,k-1
       e => e%next
     END DO
     ! second element
     ll = ns(7,i)       !element + face
     j1 = MOD(ll,10)    !face1
     k1 = (ll-j1)/10    !element1
       !  point to element k1
     IF( k1 > jelem )THEN
       n = jelem
     ELSE
       ej => heade
       n = 1
     END IF
     DO jelem=n,k1-1
       ej => ej%next
     END DO
     !  store information (j,k)  (j1,k1)
     IF( e%lnods(face(1,j)) ==  ej%lnods(face(2,j1)) .AND.  &
         e%lnods(face(2,j)) ==  ej%lnods(face(1,j1))) THEN    !lower == lower
         e%lnods(j+6)   = ej%lnods(j1) !opposite global node
         e%lnods(j+9)   = ej%lnods(j1+3) !opposite global node
         ej%lnods(j1+6) = e%lnods(j)   !opposite global node
         ej%lnods(j1+9) = e%lnods(j+3)   !opposite global node
     ELSE  !lower == upper
         e%lnods(j+6)   = ej%lnods(j1+3) !opposite global node
         e%lnods(j+9)   = ej%lnods(j1)   !opposite global node
         ej%lnods(j1+6) = e%lnods(j+3)   !opposite global node
         ej%lnods(j1+9) = e%lnods(j)     !opposite global node
     END IF
   END IF
 END DO
 DEALLOCATE (ns)
 lface = .TRUE.              !topological array created
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE toar05
