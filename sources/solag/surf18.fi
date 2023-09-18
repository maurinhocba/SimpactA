 SUBROUTINE surf18 ( heade, nelem)
 !******************************************************************
 !
 !*** get the boundary definition from the element set
 !    segments are oriented (outward normal)
 !
 !******************************************************************
 USE surf_db
 IMPLICIT NONE

 !dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nelem   !number of elements
 TYPE (ele18), POINTER :: heade          !pointer to first element in the set

 !TYPE (cont_srf), POINTER :: surfa   !INTENT(OUT) surface data

 ! local variables
 INTEGER (kind=4), PARAMETER :: nface=6 , nns=4
 INTEGER (kind=4) :: in,ie,je,jn,nf(nns),nf2(nns),ne
 INTEGER  :: nn(nns,nface), cn(nns)
 INTEGER (kind=4), ALLOCATABLE :: lface(:,:)
 TYPE (srf_seg), POINTER :: seg
 LOGICAL :: flag
 TYPE (ele18), POINTER :: e,e1

 nn = RESHAPE( (/ 1,4,3,2, 1,2,6,5, 2,3,7,6,  &
                  3,4,8,7, 4,1,5,8, 5,6,7,8 /), (/ 4,6 /) )
 cn = (/ 2,3,4,1 /)

 ALLOCATE (lface(nface,nelem))  !this arrays stores -1 if inner face
                                !  > 1 if part of the boundary
 lface = 0                      ! initializes to unknown

 !GENERATE SURFACE data_base Initialize empty list
 IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
 NULLIFY(surfa%head,surfa%tail)

 !   SEARCH FOR THE 'TWIN' FACES
 ne = 0                     !initializes number of boundary faces
 e => heade
 DO ie=1,nelem              !for each element
   jsearch : DO in=1,nface        !for each face in the element
     IF( lface(in,ie) /=0 )CYCLE  !if processed face then cycle
     nf=e%lnods(nn(:,in))        !nodes of the face
     !   search for the 'twin' face defined by nf(1:nns)
     e1 => e%next
     DO je=ie+1,nelem  !for each remaining element
       DO jn=1,nface                 !for each face
         IF( lface(jn,je) /= 0 )CYCLE  !if processed face then cycle
         nf2=e1%lnods(nn(:,jn))        !nodes of the face
         CALL compare_faces(nns,nf,nf2,flag)
         IF( flag )THEN
           lface(in,ie) = -1         !set to Inner face
           lface(jn,je) = -1         !set to Inner face
           CYCLE jsearch             !TWIN face found, exit search
         END IF
       END DO
       e1 => e1%next
     END DO
     ne = ne + 1         !increase number of boundary faces
     lface(in,ie) = ne   !assign number to boundary face
     ! store segment in list
     CALL new_seg(seg)                   !allocate surface segment
     seg%nodes(1:nns) = nf(1:nns)        !
     !WRITE(lures, '(10x,4i8)',ERR=9999) seg%nodes(1:nns)
     CALL add_seg(seg, surfa%head, surfa%tail)  !add segment to the d_b
   END DO jsearch
   e => e%next
 END DO

 surfa%nelem = ne                        !keep number of segments
 DEALLOCATE(lface)                      !release memory

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE surf18
