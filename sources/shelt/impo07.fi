 SUBROUTINE impo07 (nelem, heade, taile, ngaus, nnode, stype, ngamm, nstre)
 !******************************************************************
 !
 !*** READ element DATA for shell elements BST (ULF)
 !
 !******************************************************************
 USE meshmo_db, ONLY : nodset,numpo
 USE gvar_db, ONLY : fimpo, lab1, renum, seque, inter, overw
 IMPLICIT NONE

 INTEGER (kind=4) :: nelem,ngaus,nnode,stype,ngamm,nstre
 TYPE(ele07),POINTER:: heade,taile

 LOGICAL  found !,dummy
 INTEGER (kind=4) :: i,j,k,l,iel,chnode,isec,nlayr,osec,jsec,ln(nnode)
 REAL (kind=8) :: angle
 TYPE (ele07), POINTER :: e
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
     CALL new_ele07e (e)           !reserve memory for data
     ALLOCATE( e%dvolu(ngaus),e%lnods(nnode),e%jacin(2,2,ngaus), &
               e%stra0(6,ngaus),e%gamm0(ngamm),e%strsg(nstre,ngaus),     &
               e%qbar(ngamm),e%ambda(2,ngaus))
     ALLOCATE (e%lnods(nnode))
     IF(ngaus == 3)THEN
       ALLOCATE(e%cartd(nnode,2,ngaus))
     ELSE
       ALLOCATE(e%cartd(3,2,ngaus))
       IF(nnode == 9) ALLOCATE(e%cd(4,2,3))
     END IF
     READ (fimpo) e%numel,e%matno,e%lnods(1:nnode),e%angle
     ! modify nodes
     IF (renum)THEN
       IF( seque )THEN
         DO j=1,nnode
           k = e%lnods(j)
           IF( k == 0 )CYCLE
           DO l=1,numpo
             IF( k == nodset(l) )THEN
               e%lnods(j) = lab1 + l
               EXIT
             END IF
           END DO
         END DO
       ELSE
         DO j=1,nnode
           k = e%lnods(j)
           IF( k == 0 )CYCLE
           e%lnods(j) = e%lnods(j) + lab1
         END DO
       END IF
     END IF
     !     original nodes are changed
     DO i=1,6
       e%lnods(i) = chnode(e%lnods(i))
     END DO
     DO i=7,nnode
       IF(e%lnods(i) > 0) e%lnods(i) = chnode(e%lnods(i))
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
       CALL sect_search (jsec,found,sec,k)     !search in list of sections
       jsec = k
     END IF
     e%matno = jsec           !associated section
     IF( nstre > 8 .AND. sec%secty == 15) ALLOCATE(e%nangl(2,nnode),e%jas(2,ngamm))

   END IF
   nlayr = sec%iprop(1)         !number of layers (integration points)
  ! get internal variables
   IF(inter) THEN
     READ(fimpo) l
     IF( .NOT. overw )ALLOCATE( e%ehist(5,ngaus))
     IF( l == 1 )THEN
       READ(fimpo) ((e%ehist(i,j),i=1,5),j=1,ngaus)
       IF( stype == 0 )THEN  !layered large strain plasticity
           IF( .NOT. overw )ALLOCATE( e%strap(6*nlayr,ngaus),e%stres(5,nlayr,ngaus) )
           READ(fimpo) ((e%strap(i,j),i=1,6*nlayr),j=1,ngaus)
           e%stres = 0d0
         ELSE   !solid isotropic
           IF( .NOT. overw )ALLOCATE( e%strap(nstre,ngaus))
           READ(fimpo) ((e%strap(i,j),i=1,nstre),j=1,ngaus)
           NULLIFY(e%stres)
         END IF
     ELSE
       NULLIFY(e%ehist,e%strap)
     END IF
   END IF
   e%qbar = 0d0
   e%strsg = 0d0
   IF( overw )THEN
     e => e%next
   ELSE
     CALL add_ele07e(e, heade, taile)
   END IF
 END DO

 RETURN

 END SUBROUTINE impo07
