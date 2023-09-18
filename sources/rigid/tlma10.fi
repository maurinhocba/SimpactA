 SUBROUTINE tlma10 (nelem,ntype,ndime,nnode,ngaus,lnods,dvolu,     &
                    shape,matno,ymass,iftmp)
 !
 ! computes nodal mass for a set of rigid elements
 ! compute capacity-mass vector for element RIGID solid element
 !
 USE npo_db, ONLY : tempe
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem, & !number of elements
                                 ntype, & !type of rigid body
                                 ndime, & !problem dimension
                                 nnode, & !number of nodes per element
                                 ngaus, & !number of gauss points
                                 lnods(:,:), &  !connectivities
                                 iftmp(:,:)     !associated dofs
 INTEGER (kind=4), POINTER :: matno(:)      !properties
 REAL (kind=8), INTENT (IN OUT) :: ymass(:), & !nodal mass
                                   dvolu(:,:), & !Gauss point volume
                                   shape(:,:)    !shape functions
 INTEGER (kind=4) :: i,j,k,g

 TYPE(section), POINTER :: sec
 TYPE (curve), POINTER :: c
 INTEGER (kind=4) isec,iel,osec,numpt
 REAL (kind=8) :: rhoel,thick,voli,heatc,val,fac
 REAL (kind=8) :: diagm(nnode)

 osec = -1
 DO iel=1,nelem
   isec = matno(iel)                    !associated material
   IF( isec /= osec )THEN
     sec => psecs(isec)%p
     rhoel = sec%mtbas%prope(5)          !material density
     heatc = sec%mtbas%prope(7)          !material capacity
     numpt = INT(sec%mtbas%props(7))     !associated curve for Heat capacity
     CALL cur_point( numpt,sec%mtbas%chead,c) !point to the curve
     thick = 1d0
     IF(ndime == 2 .AND. ntype == 1 )thick = sec%rprop(1)      ! thickness
     fac = thick*rhoel*heatc
     osec = isec
   END IF
   IF( numpt > 0 )THEN !compute Heat capacity as a function of temperature
     ! temperature at element center (middle surface)
     val = SUM(tempe(1,lnods(1:nnode,iel)))/nnode
     i = 1
     heatc = inte_cr(c%val,c%np,val,i)
     fac = thick*rhoel*heatc
   END IF
   ! ---computes jacobian at Gauss point
   diagm = 0d0
   DO g=1,ngaus
     voli=dvolu(g,iel)*fac  !element capacity
     ! ---generates lumped matrix proportional to diagonal
     DO j=1,nnode !for each node
       diagm(j) = diagm(j) + voli*shape(j,g)
     END DO
   END DO
   ! ---add to global CAPACITY
   DO i=1,nnode                          !for each element node
     j=lnods(i,iel)                  !node number
     k=iftmp(1,j)                !DOF
     IF( k > 0 )ymass(k) = ymass(k) + diagm(i)
   END DO

 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE tlma10
