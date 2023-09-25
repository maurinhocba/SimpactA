 SUBROUTINE surf19 ( heade, nelem)
 !******************************************************************
 !
 !*** get the boundary definition from the element set
 !    segments are oriented (outward normal)
 !
 !******************************************************************
 USE surf_db , ONLY : srf_seg,delete_seg,add_seg,new_seg,surfa
 IMPLICIT NONE
 !dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nelem      !number of elements
 TYPE (ele19), POINTER :: heade   !pointer to first element in the set
 !TYPE (cont_srf), POINTER :: surfa   !INTENT(OUT) surface data

 ! local variables
 INTEGER (kind=4) :: ii,ie,n1,n2,je,jn,m1,m2,ne
 INTEGER (kind=4), ALLOCATABLE :: lside(:,:),order(:,:)
 TYPE (srf_seg), POINTER :: seg
 TYPE (ele19), POINTER :: e,e1

 ALLOCATE (lside(nnode,nelem))  !this arrays stores -1 if inner side
                                !  > 1 if part of the boundary
 lside = 0                      ! initializes to unknown

 !   SEARCH FOR THE 'TWIN' SIDES
 ne = 0                     !initializes number of boundary sides
 e => heade
 DO ie=1,nelem              !for each element
   jsearch : DO ii=1,nnode        !for each side in the element
     IF( lside(ii,ie) /=0 )CYCLE  !if processed side then cycle
     n1=e%lnods(ii)          !first node of the side
     n2=e%lnods(MOD(ii,nnode)+1)      !second node of the side
     IF( n1 == n2 ) THEN      !if both node are the same ==> not a side
       lside(ii,ie) = -1      !set as if where an inner side
       CYCLE                  !next side
     END IF
     !   search for the 'twin' side of n1-n2
     e1 => e%next
     DO je=ie+1,nelem  !for each remaining element
       DO jn=1,nnode                 !for each side
         IF( lside(jn,je) /= 0 )CYCLE   !if processed side then cycle
         m1=e1%lnods(jn)                  !first node of the side
         m2=e1%lnods(MOD(jn,nnode)+1)     !second node of the side
         IF( m1 == n2 .AND. m2 == n1) THEN  !if this is the TWIN side
           lside(ii,ie) = -1         !set to Inner side
           lside(jn,je) = -1         !set to Inner side
           CYCLE jsearch             !TWIN side found, exit search
         END IF
       END DO
       e1 => e1%next
     END DO
     ne = ne + 1         !increase number of boundary sides
     lside(ii,ie) = ne   !assign number to boundary side
   END DO jsearch
   e => e%next
 END DO

 ! SORT sides according to connectivities
 ALLOCATE ( order(3,ne) )
 order = 0   !initializes
 e => heade
 DO ie=1,nelem              !for each element
   DO ii=1,nnode                !for each side in the element
     IF( lside(ii,ie) < 0 )CYCLE  !if an inner side then cycle
     n1=lside(ii,ie)          !segment order
     n2=e%lnods(ii)           !second node of the side
     !   search for the 'next' side of n1-n2
     e1 => heade
     nsearch : DO je=1,nelem  !for each remaining element (including IE)
       DO jn=1,nnode                 !for each side
         IF( lside(jn,je) <  0 )CYCLE   !if an inner side then cycle
         m1=e1%lnods(MOD(jn,nnode)+1)   !first node of the side
         IF( n2 == m1) THEN  !if this is the next side
           m2=lside(jn,je)            !segment order
           order(1,n1) = ii           !node   (this side)
           order(2,n1) = ie           !element(this side
           order(3,n1) = m2           !node   (next side in order)
           EXIT nsearch               !next segment found, exit search
         END IF
       END DO
       e1 => e1%next
     END DO nsearch
   END DO
   e => e%next
 END DO

 !GENERATE SURFACE data_base
 !Initialize empty list
 IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
 NULLIFY(surfa%head,surfa%tail)
 ii = 1
 DO ie=1,ne     !for each side
   IF(order(1,ii) == 0)THEN  !find a non-stored side
     ii = 1         !initializes
     DO
       IF(order(1,ii) /= 0)EXIT  !non-stored side
       ii = ii+1                 !next side
     END DO
   END IF
   jn = order(1,ii)     !side order in element connectivities
   je = order(2,ii)     !element
   e => heade
   DO n1=1,je-1
     e => e%next
   END DO
   ! store segment in list
   CALL new_seg(seg)                         !allocate surface segment
   seg%nodes(2) = e%lnods(jn)                !first node
   seg%nodes(1) = e%lnods(MOD(jn,nnode)+1)   !second node
   !WRITE(lures, '(10x,2i8)',ERR=9999) seg%nodes(1:2)
   CALL add_seg(seg, surfa%head, surfa%tail)  !add segment to the d_b
   order(1,ii) = 0      !set as stored or non-existent
   ii = order(3,ii)     !next side
 END DO
 surfa%nelem = ne                       !keep number of segments
 DEALLOCATE(lside, order)               !release memory

 RETURN
 END SUBROUTINE surf19
