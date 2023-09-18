 SUBROUTINE surf10 (lnods,nelem,nnode,ndime,ntype)
 !
 !  Gets surface for contact
 !  output trought SURF_DB
 !
 USE surf_db
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: lnods(:,:), & ! connectivities
                                 nelem,      & ! number of elements
                                 nnode,      & ! number of nodes per element
                                 ntype,      & ! rigid body type
                                 ndime         ! problem dimensnion

 INTEGER (kind=4) :: ie,i
 TYPE (srf_seg), POINTER :: seg
 ! connectivities from 6-node triangles to 3-node triangles (NT)
 !            and from 3-node lines to 2-node lines
 INTEGER (kind=4), PARAMETER :: nt(3,4) = (/ 1,4,6, 2,5,4, 3,6,5, 4,5,6 /), &
                                ns(2,2) = (/ 1,2, 2,3 /)

 INTERFACE
   INCLUDE 'surf03.h'
   INCLUDE 'surf05.h'
 END INTERFACE

 SELECT CASE (ntype)
 CASE(0)   !rigid body defined by particles
   IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
   NULLIFY(surfa%head,surfa%tail)          !Initialize empty list
   !CALL new_seg(seg)                         !allocate surface segment
   !seg%nodes(1) = mastn                      !master node
   !CALL add_seg(seg, surfa%head, surfa%tail) !add segment to the d_b
   DO ie=1,nnode     !for each node
     ! store nodes in list
     CALL new_seg(seg)                          !allocate surface segment
     seg%nodes(:) = (/ lnods(ie,1), 0, 0, 0 /) !node
     CALL add_seg(seg, surfa%head, surfa%tail) !add segment to the d_b
   END DO
   surfa%nelem = nnode !+ 1                    !keep number of nodes
 CASE(1:3) !solid rigid body
   IF(ndime == 2)THEN
     IF( nnode <= 4 )THEN !bilinear quads or linear triangles
       CALL surf03 (lnods,nelem,nnode)
     ELSE
       CALL runen3('SURF10: cannot obtain 8/9-node elements boundary')
     END IF
   ELSE  !IF (ndime == 3) trilinear hexahedra or linear tetrahedra
     CALL surf05 (lnods,nelem,nnode)
   END IF
 CASE(4:6) !rigid surface
   !GENERATE SURFACE data_base
   IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
   NULLIFY(surfa%head,surfa%tail)          !Initialize empty list
   DO ie=1,nelem     !for each element
     ! store segment in list
     IF( ndime ==2 )THEN                        !2-d problems
       SELECT CASE (nnode)
       CASE (2)                                 !linear elements
         CALL new_seg(seg)                          !allocate surface segment
         seg%nodes(1:nnode) = lnods(:,ie)           !nodes
         !WRITE(lures, '(10x,6i8)') seg%nodes(1:nnode)
         CALL add_seg(seg, surfa%head, surfa%tail) !add segment to the d_b
       CASE (3)                                 !cuadratic elements
         DO i=1,2                               !for each 2-node segment
           CALL new_seg(seg)                             !allocate surface segment
           seg%nodes = (/ lnods(ns(:,i),ie), 0, 0 /)  !nodes
           !WRITE(lures, '(10x,6i8)') seg%nodes(1:2)
           CALL add_seg(seg, surfa%head, surfa%tail)  !add segment to the d_b
         END DO
       END SELECT
     ELSE                                       !3-d problems
       SELECT CASE (nnode)
       CASE (3,4)                               !linear elements (quads and triangles)
         CALL new_seg(seg)                             !allocate surface segment
         seg%nodes(1:nnode) = lnods(:,ie)           !nodes
         !WRITE(lures, '(10x,6i8)') seg%nodes(1:nnode)
         CALL add_seg(seg, surfa%head, surfa%tail) !add segment to the d_b
       CASE (6)                                 !6-node triangle
         DO i=1,4                               !for each linear triangle
           CALL new_seg(seg)                           !allocate surface segment
           seg%nodes = (/ lnods(nt(:,i),ie), 0 /)           !nodes
           !WRITE(lures, '(10x,6i8)') seg%nodes(1:3)
           CALL add_seg(seg, surfa%head, surfa%tail)  !add segment to the d_b
         END DO
       END SELECT
     END IF
   END DO
   IF( nnode == 6 ) THEN
     surfa%nelem = 4*nelem                   !keep number of segments
   ELSE IF(ndime == 2 .AND. nnode == 3 ) THEN
     surfa%nelem = 2*nelem                   !keep number of segments
   ELSE
     surfa%nelem = nelem                     !keep number of segments
   END IF
 END SELECT

 RETURN
 END SUBROUTINE surf10
