 SUBROUTINE impo18 (nelem, ngaus, heade, taile, shell)
 !******************************************************************
 !
 !*** READ element DATA for Solid elements SOLAG (TLF)
 !
 !******************************************************************
 USE meshmo_db, ONLY : nodset,numpo
 USE gvar_db, ONLY : fimpo, lab1, renum, seque, inter, overw
 IMPLICIT NONE

 INTEGER (kind=4), PARAMETER ::  nnode=8
 INTEGER (kind=4) :: nelem,ngaus
 TYPE(ele18),POINTER:: heade,taile
 LOGICAL :: shell

 LOGICAL  found !,dummy
 INTEGER (kind=4) :: i,j,k,l,iel,chnode,isec,osec,nvare,jsec,nv,ln(nnode)
 REAL (kind=8) :: angle(3)
 TYPE (ele18), POINTER :: e
 TYPE(section), POINTER :: sec

 !
 ! *** READ the element nodal connections, and the property numbers.
 !
 osec = -1
 iel = 0
 IF( overw ) e => heade
 DO iel=1,nelem
   ! loop over elements
   IF( overw )THEN
     READ (fimpo)i,j,ln,angle
     IF( i /= e%numel )WRITE(55,"(' import set element label changed from ', &
          &                         i6,' to ',i6)")e%numel,i
     !  section must be modified
     isec = e%matno                 !associated material
     IF( isec /= osec )  sec => psecs(isec)%p  !point to new section
   ELSE
     CALL new_ele18e(e)           !reserve memory for data
     ALLOCATE (e%cartd(nnode,3,ngaus), e%dvol(ngaus), e%stint(6,ngaus))    !reserve memory for data
     READ (fimpo) e%numel,e%matno,e%lnods(1:nnode),e%angle
     IF( shell )THEN
       ALLOCATE(e%nfdas(8,4),e%jacin(2,2,ngaus))
     ELSE
       NULLIFY( e%nfdas,e%jacin )
     END IF

     ! modify nodes
     IF (renum)THEN
       IF( seque )THEN
         DO j=1,nnode
           k = e%lnods(j)
           DO l=1,numpo
             IF( k == nodset(l) )THEN
               e%lnods(j) = lab1 + l
               EXIT
             END IF
           END DO
         END DO
       ELSE
         e%lnods(1:nnode) = e%lnods(1:nnode) + lab1
       END IF
     END IF
     !     original nodes are changed
     DO i=1,nnode
       e%lnods(i) = chnode(e%lnods(i))
     END DO
     !  section must be modified
     isec = e%matno                 !associated material
     IF( isec /= osec )THEN         !if a different section than previous
       jsec = 1
       DO
         IF(snn(1,jsec) == isec )EXIT
         jsec = jsec+1
       END DO
       jsec = snn(2,jsec)           !associated section
       jsec = snn(2,jsec)           !associated section
       CALL sect_search (jsec,found,sec,k)     !search in list of sections
       jsec = k
     END IF
     e%matno = jsec           !associated section

   END IF
   IF( isec /= osec )THEN         !if a different section than previous
     osec = isec                  !keep section label
     nvare = sec%iprop(2)         !number of internal variables per layer
   END IF
  ! get internal variables
   IF(inter) THEN
     READ(fimpo) nv
     IF( nv > 0 ) THEN
       IF( nv > nvare ) STOP 'error in NVARE: number of internal variables'
       IF( .NOT. ASSOCIATED(e%gausv)) ALLOCATE (e%gausv(nvare,ngaus))
       READ(fimpo) ((e%gausv(i,j),i=1,nv),j=1,ngaus)
     ELSE
       IF(ASSOCIATED(e%gausv)) DEALLOCATE (e%gausv)
     END IF
   END IF
   IF( overw )THEN
     e => e%next
   ELSE
     CALL add_ele18e(e, heade, taile)
   END IF
 END DO

 RETURN

 END SUBROUTINE impo18
