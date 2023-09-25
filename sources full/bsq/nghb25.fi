 SUBROUTINE nghb25(lside,nelem,nbs,heade,bhead)

 ! Create topological array

 IMPLICIT NONE
 LOGICAL :: lside
 INTEGER(kind=4),  INTENT(IN):: nelem  !number of elements
 INTEGER (kind=4), INTENT(OUT) :: nbs  !number of branching sides
 TYPE (ele25), POINTER :: heade   !pointer to first element in the set
 TYPE (sideb), POINTER :: bhead   !pointer to first branching side

 ! local variables
  TYPE (ele25), POINTER  :: e, ej                       !element pointers

 INTEGER (kind=4) maxsd,i,ii,j,jj,k,kk,l,n,nsd,j1,k1,ielem,jelem,nn,jjj
 !INTEGER (kind=4), PARAMETER :: nextn(3) = (/2,3,1/) !cycling list
 LOGICAL :: first
 TYPE (sideb), POINTER :: si,btail  !pointer to a branching side
 INTEGER (kind=4), ALLOCATABLE :: ns(:,:)

 !_____________________________
 maxsd = 3*nelem+1    !maximum number of sides
 ALLOCATE (ns(20,maxsd))   !auxiliar array (addmits 8 elements over one side)
 ns = 0                 !initializes array
 nsd = 0                !initializes number of sides

 e => heade             !point to first element
  DO i=1,nelem          !for each element
    DO ii=1,4           !for each side in the element
      j = e%lnods(ln(2,ii))     !internal node numeration (first node of the side)
      k = e%lnods(ln(1,ii))     !internal node numeration (second node)
      first = j < k             !J is the node with lower numeration
      l = MIN(j,k)              !use the lower value
      k = MAX(j,k)              !greater value
      j = l                     !lower value
      jj = 1                    !initializes pointer
      DO                  !search in existing sides
        IF( jj > nsd )THEN  !all sides checked ?
          nsd = jj          !increase number of sides
          ns(1,jj) = j      !keep side nodes
          ns(2,jj) = k
          ns(3,jj) = 4      !position
          EXIT
        ELSE IF(ns(1,jj) == j .AND. ns(2,jj) == k )THEN  !side exists
          ns(3,jj) = ns(3,jj) + 1   !position
          EXIT
        END IF
        jj = jj + 1         !next side
      END DO
     ! keep element and side in one value
     l = ns(3,jj)          !position in array NS
     IF( first )THEN   !ordered
       ns(l,jj) =   ii+10*i    !keeps ii = side in element and i = element
     ELSE              !reverse order
       ns(l,jj) = -(ii+10*i)
     END IF
   END DO
   e => e%next          !go to next element
 END DO
 ! generates data base for branching nodes
 e => heade             !point to first element
 ielem = 1              !initializes element I
 ej => heade            !point to first element
 jelem = 1              !initializes element J
 nbs = 0                !initializes number of branching nodes
 DO i=1,nsd             !for each side
   nn = ns(3,i)-3       !number of connected elements
   IF( nn > 2) THEN   !if a branching side
     nbs = nbs + 1        !increase number
     ALLOCATE (si)        !generates pointer
     IF(nbs == 1)THEN   !for first branching side
       bhead => si      !keep first pointer
       btail => si      !auxiliar last pointer
     ELSE
       btail%next => si !updates pointer to new branching side
       btail => si      !updates pointer
     END IF
     si%nn = nn         !keep number of connected elements
     ALLOCATE( si%lnods(-1:2*nn), & ! j,k, i=1,2*nn (-) inverse order in conns
               si%alph0(nn),    & ! i=1,nn-1 angles between elems. nn = 1/ls
               si%fc(nn,nn),    & ! (1:nn,1:nn-1) matrix F  (1:nn,nn) relative stiffness
               si%t0(3,nn),  &
               si%tb(3,nn),  &
               si%c(4,nn,nn),      & ! normal derivative of shape functions
               si%gamma(nn),    & ! distortion at each element
               si%bbr(3,-1:2*nn) )   ! average rotation B - matrix

     si%lnods(-1:0) = ns(1:2,i)    !side nodes j-k (ordered)
     kk = 0
     DO ii=1,nn                     !for each element on this side
       !kk = kk + ii
       jj = ns(ii+3,i)     !element + side
       j = MOD(ABS(jj),10) !side
       k = ABS(jj-j)/10    !element
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
       !  transfer information
       e%si(j)%p => si     !keep a pointer to the branching side
       !e%lnods(j+4) = -ii   !keep position in side connectivities
       e%lnods(ln(3,j)) = -ii
       e%lnods(ln(4,j)) = -ii
       IF( jj > 0 )THEN              !rigth direction
         si%lnods(ii+kk) =  e%lnods(nextn(j+1))  !keep node & order (+sign)
         si%lnods(ii+kk+1) =  e%lnods(nextn(j+2))
       ELSE
         si%lnods(ii+kk) = -e%lnods(nextn(j+1))  !keep node & order (-sign)
         si%lnods(ii+kk+1) = -e%lnods(nextn(j+2))
       END IF
       e%lside(j) = -nbs   !unnecessary perhaps
       si%fc(ii,nn) = REAL(ielem)    !keep element number
       kk = kk + 1
     END DO

   ELSE IF( nn == 2 )THEN          !standard element

     !  first element
     jj = ns(4,i)          !element + side
     j = MOD(ABS(jj),10)   !side
     k = ABS(jj-j)/10      !element
     jjj= jj/k           !+/-1
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
     kk = ns(5,i)          !element + side
     j1 = MOD(ABS(kk),10)  !side1
     k1 = ABS(kk-j1)/10    !element1
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
    !  store information
    ! IF( jjj*kk > 0 )THEN
    !   WRITE(55,"('  WARNING, non-orientated surface between '/, &
    !              & '  elements ',i6,' and ',i6)",ERR=9999) e%numel,ej%numel
    ! END IF
     e%lnods(fn(j))    = ej%lnods(nextn(j1+1)) ! global node
     e%lnods(fn(j)+1)  = ej%lnods(nextn(j1+2)) ! global node
     e%lside(j)        = jelem                    ! adjacent element
     ej%lnods(fn(j1))  = e%lnods(nextn(j+1))  ! global node
     ej%lnods(fn(j1)+1)= e%lnods(nextn(j+2))  ! global node
     ej%lside(j1)      = ielem                    ! adjacent el
   END IF
 END DO
 DEALLOCATE (ns)
 lside = .TRUE.              !topological array created

 RETURN
 END SUBROUTINE nghb25
