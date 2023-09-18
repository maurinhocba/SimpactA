 SUBROUTINE toar11(heade,lside,nbn,nhead,nelem)

 ! Create topological array

 USE ctrl_db, ONLY : npoin
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem
 INTEGER (kind=4), INTENT(OUT) :: nbn
 TYPE (ele11), POINTER :: heade   !pointer to first element in the set
 TYPE (nodei), POINTER :: nhead   !pointer to first branching node
 LOGICAL, INTENT(OUT) :: lside

 ! local variables
 TYPE (ele11), POINTER  :: e, e1    !element pointers
 INTEGER (kind=4) i,j,k,l,m,nn,osec,isec,jsec
 REAL (kind=8) :: ei,ej
 TYPE (nodei), POINTER :: ni,ntail  !pointer to a branching node
 TYPE (section), POINTER :: sec,secj
 LOGICAL :: found
 INTEGER (kind=4), ALLOCATABLE :: np(:)

 !_____________________________
 ALLOCATE (np(npoin))   !auxiliar array
 np = 0                 !initializes
 e => heade             !point to first element
 ! count appearence of each node
 DO i=1,nelem           !for each element
   DO j=2,3             !2-node connectivities
     k = e%lnods(j)     !internal node numeration
     np(k) = np(k) + 1  !increase number of appearences
   END DO
   e%lnods(1)= 0        !set to 0
   e%lnods(4)= 0        !set to 0
   e => e%next          !go to next element
 END DO
 ! generates data base for branching nodes
 nbn = 0                !initializes number of branching nodes
 DO i=1,npoin           !for each node
   IF(np(i) > 2) THEN   !if a branching node
     nbn = nbn + 1      !increase number
     ALLOCATE (ni)      !generates pointer
     IF(nbn == 1)THEN   !for first branching node
       nhead => ni      !keep first pointer
       ntail => ni      !auxiliar last pointer
     ELSE
       ntail%next => ni !updates pointer to new branching node
       ntail => ni      !updates pointer
     END IF
     nn = np(i)         !number of connected elements
     ni%nn = nn         !keep this value
     ALLOCATE( ni%lnods(0:nn), ni%alph0(nn), ni%fc(nn), ni%bb(2,0:nn) ) !allocates memory
     ni%lnods(0) = i
     ni%lnods(1:nn) = 0
     np(i) = -nbn       !keep position in list
   END IF
 END DO
 !_____________________________

 osec = -1
 e => heade             !point to first element
 DO i=1,nelem
   isec = e%matno
   IF( isec /= osec )THEN
     sec => psecs(isec)%p
     ei = sec%rprop(12)
     osec = isec
   END IF
   DO j=2,3                             !for each node in the element
     k = e%lnods(j)                     !node number
     SELECT CASE ( np(k) )              !node type
     CASE (:-1)            !branching node
       ni => nhead         !point to top of data base
       DO l=1,-np(k)-1     !loop to the correct position
         ni => ni%next
       END DO
       l = 1               !initializes
       DO                  !loop until a free value
         IF(ni%lnods(l) == 0 )THEN
           ni%lnods(l) = e%lnods(nyn(j)) !include opposite node in the data base
           e%lnods(nxn(j)) = np(k)
           ni%fc(l) = ei
           EXIT                          !exit loop
         END IF
         l = l+1
       END DO
       e%lnods(j+3) = -2  !branching node
     !CASE (0)    !not possible
     !CASE (1)    !end node, do not search
     !  CYCLE
     CASE (2)    !search for node
       IF( e%lnods(nxn(j)) /= 0 )CYCLE !node already find
       e1 => e%next
       search : DO l=i+1,nelem
         DO m=2,3                             !for each node in the element
           IF( e1%lnods(m) == k )THEN         !nodes are coincidents
             e%lnods(nxn(j)) = e1%lnods(nyn(m))   !
             e%lnods(j+3) = l                     !keep neighbor element
             e1%lnods(nxn(m)) = e%lnods(nyn(j))
             e1%lnods(m+3) = i                    !keep neighbor element
             jsec = e1%matno
             IF( isec /= jsec )THEN
               secj => psecs(jsec)%p
               ej = secj%rprop(12)
             ELSE
               ej = ei
             END IF
             IF( j == 2 )THEN
               e%cab(4,1) = ej
             ELSE
               e%cab(4,2) = ej
             END IF
             IF( m == 2 )THEN
               e1%cab(4,1) = ei
             ELSE
               e1%cab(4,2) = ei
             END IF
             EXIT search
           END IF
         END DO
         e1 => e1%next
       END DO search
     END SELECT
   END DO
   e => e%next               !point to next element to process
 END DO
 lside = .TRUE.              !topological array created
 RETURN
 END SUBROUTINE toar11
