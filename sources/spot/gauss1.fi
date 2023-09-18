 SUBROUTINE gauss1(ndime,nnode,neulr,nelem,heade,coord,eule0,lowers,uppers,gauss)
 !********************************************************************
 !
 !*** evaluates gauss points constants for truss element
 !
 !********************************************************************
 USE npo_db, ONLY : label
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: ndime,nnode,neulr,nelem
 REAL    (kind=8), INTENT(IN) :: coord(:,:),eule0(:,:)
 TYPE (ele01), POINTER :: heade
 CHARACTER (len=*), INTENT(IN OUT) :: lowers,uppers
 LOGICAL, INTENT(IN OUT) :: gauss

 INTEGER (kind=4) ielem,n1,n2,nm,ln(6),i,j
 TYPE (ele01), POINTER :: e
 LOGICAL :: same,sur,sfound,efound
 REAL(kind=8) :: x(3,nnode),xc(3),t(3,3),long,a2,a2p,l1(3),l2(3),y(3),z(3), &
                 z1,z2,xa(3),l1p(3),l2p(3),rm(9)

 INTERFACE
   INCLUDE 'inrotm.h'
 END INTERFACE


 sur = LEN_TRIM(lowers) > 0
 e => heade
 DO ielem = 1,nelem
   IF( nnode == 2 )THEN
     n1 = e%lnods(1)
     n2 = e%lnods(2)
     same = (n2 == 0)
     IF( same )THEN !keep original coordinates
       e%gausv(2:1+ndime) = coord(1:ndime,n1)
       nm = 3+ndime
       IF( e%rot)THEN
         IF( neulr == 1 )THEN
           e%gausv(nm) = eule0(1,n1)
         ELSE
           CALL inrotm(eule0(:,n1),rm(:))
           e%gausv(nm:nm+neulr-1) = rm
         END IF
       END IF
     END IF
   ELSE     !
     IF( sur )THEN
       nm = e%lnods(1); n1 = 0; n2 =0  !auxiliar node
       CALL search_ele13 (ln(1:3), lowers, n1, sfound, efound, nm, coord, .TRUE.)
       IF( .NOT.sfound )THEN             !set given but not found
         WRITE (lures, "(' Element set not found ')",ERR=9)
         CALL runend('lower surface not found: invalid data ')
       ELSE IF( .NOT.efound )THEN        !element not found
         WRITE (lures, "(' Element not found in the set ')",ERR=9)
         CALL runend('no proyection on lower surface')
       END IF
       CALL search_ele13 (ln(4:6), uppers, n2, sfound, efound, nm, coord, .TRUE.)
       IF( .NOT.sfound )THEN             !set given but not found
         WRITE (lures, "(' Element set not found ')",ERR=9)
         CALL runend('upper surface not found: invalid data ')
       ELSE IF( .NOT.efound )THEN        !element not found
         WRITE (lures, "(' Element not found in the set ')",ERR=9)
         CALL runend('no proyection on upper surface')
       END IF
       e%lnods(1:6) =  ln   !keep nodes
       !write(55,"(6i5)")label(ln)
     ELSE
       nm = 0
       IF( e%numel < 0 ) THEN
         nm = -e%numel
         e%numel = ielem !modify label
       END IF
     END IF
     x(:,1:6) = coord(:,e%lnods(1:6))                     !coordinates
     IF(nm == 0)THEN
       xa = (x(:,1)+x(:,2)+x(:,3))/3d0     !no auxiliar node
     ELSE
       xa = coord(:,nm)
     END IF
     ! compute local system on lower element
     l1 = x(:,1)-x(:,3)                       !second side
     l2 = x(:,2)-x(:,3)                       !first side (-)
     CALL vecpro(l1,l2,t(1,3))                !normal
     t(:,1) = l1
     CALL vecuni(3,t(1,1),a2)                 !direction 1
     CALL vecuni(3,t(1,3),a2)                 !unit normal and twice the area
     CALL vecpro(t(1,3),t(1,1),t(1,2))        !direction 2
     ! compute lower element proyection
     y = xa - x(:,3)
     long = DOT_PRODUCT(y,t(:,3))
     y = y - long*t(:,3)
     CALL vecpro(y,l2,z)                      !auxiliar
     e%gausv(1) = DOT_PRODUCT(z,t(:,3))/a2    !first coordinate
     CALL vecpro(l1,y,z)                      !auxiliar
     e%gausv(2) = DOT_PRODUCT(z,t(:,3))/a2    !second coordinate
     ! compute upper element proyection
     l1p = x(:,4)-x(:,6)                       !second side
     l2p = x(:,5)-x(:,6)                       !first side (-)
     CALL vecpro(l1p,l2p,z)                     !normal
     a2p = DOT_PRODUCT(z,t(:,3))              !twice the area
     y = xa - x(:,6)                      !distance vector
     long = DOT_PRODUCT(y,t(:,3))             !normal proyection
     y = y - long*t(:,3)                      !in plane vector
     CALL vecpro(y,l2p,z)                      !auxiliar
     e%gausv(3) = DOT_PRODUCT(z,t(:,3))/a2p   !first coordinate
     CALL vecpro(l1p,y,z)                      !auxiliar
     e%gausv(4) = DOT_PRODUCT(z,t(:,3))/a2p   !second coordinate
     z1 = 1d0 - e%gausv(1) - e%gausv(2)
     z2 = 1d0 - e%gausv(3) - e%gausv(4)
     xc = -e%gausv(1)*x(:,1) -e%gausv(2)*x(:,2) - z1*x(:,3) &   !distance vector between centers
          +e%gausv(3)*x(:,4) +e%gausv(4)*x(:,5) + z2*x(:,6)
     ! compute reference values of the fastener
     e%gausv(11) = DOT_PRODUCT(xc,t(:,3))       !distance along normal
     e%gausv(12) = DOT_PRODUCT(xc,t(:,1))       !initial distortion dir 1
     e%gausv(13) = DOT_PRODUCT(xc,t(:,2))       !initial distortion dir 2
     ! compute upper node relative coordinates
     j = 5
     DO i=4,6
       y = x(:,i) - x(:,3)
       long = DOT_PRODUCT(y,t(:,3))
       y = y - long*t(:,3)
       CALL vecpro(y,l2,z)                      !auxiliar
       e%gausv(j) = DOT_PRODUCT(z,t(:,3))/a2    !first coordinate
       CALL vecpro(l1,y,z)                      !auxiliar
       e%gausv(j+1) = DOT_PRODUCT(z,t(:,3))/a2  !second coordinate
       ! compute reference value at the node
       z1 = 1d0 - e%gausv(j) - e%gausv(j+1)
       z = x(:,i) - e%gausv(j)*x(:,1) - e%gausv(j+1)*x(:,2) - z1*x(:,3)
       e%gausv(i+10) = DOT_PRODUCT(z,t(:,3)) !initial distance nod i
       j = j+2
     END DO
   END IF
   e => e%next
 END DO
 gauss = .FALSE.
 IF( sur )THEN !forget names
   lowers = ''
   uppers = ''
 END IF
 RETURN
 9 CALL runen2('')

 END SUBROUTINE gauss1
