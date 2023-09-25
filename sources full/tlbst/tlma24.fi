 SUBROUTINE tlma24 (heade,ymass,iftmp)

 ! compute capacity-mass vector for element BST-CST (TLF) shell element

 IMPLICIT NONE
 !  dummy arguments

 REAL (kind=8) :: ymass(:)        !lumped nodal mass
 INTEGER (kind=4) :: iftmp(:,:)
 TYPE (ele24), POINTER  :: heade   !pointer to first element

 ! local variables
 INTEGER (kind=4) isec,lnode,idofn,j,i,osec, &
                  numpt    !number of points in curve
 REAL (kind=8) :: diagm,voli,rhoel,thick,heatc,val
 TYPE (ele24), POINTER :: elem
 TYPE (section), POINTER :: sec
 TYPE (curve), POINTER :: c

 osec = -1
 elem => heade    !point to first element
 DO
   IF (.NOT.ASSOCIATED (elem) ) EXIT     !if all elements processed, exit
   isec = elem%matno                    !associated material
   IF( isec /= osec )THEN
     sec => psecs(isec)%p
     rhoel = sec%mtbas%prope(5)          !material density
     heatc = sec%mtbas%prope(7)          !material capacity
     numpt = INT(sec%mtbas%props(7))     !associated curve for Heat capacity
     CALL cur_point( numpt,sec%mtbas%chead,c) !point to the curve
     thick = sec%rprop(1)                !shell thickness
     osec = isec
   END IF
   IF( numpt > 0 )THEN !compute Heat capacity as a function of temperature
     ! temperature at element center (middle surface)
     val = SUM(tempe(1,elem%lnods(1:3)))/3d0
     IF( ndoft == 2 ) val = val/2d0+ SUM(tempe(2,elem%lnods(1:3)))/3d0
     i = 1
     heatc = inte_cr(c%val,c%np,val,i)
   END IF
   voli=elem%area1*thick*rhoel*heatc     !element capacity
   diagm = voli/3.d0                     !nodal mass
   !
   ! ---generates lumped mass matrix proportional to diagonal
   !
   SELECT CASE (ndoft)
   CASE (1)
     DO i=1,3                          !for each element node
       lnode=elem%lnods(i)             !node number
       idofn=iftmp(1,lnode)                !middle surface
       IF( idofn > 0 )ymass(idofn) = ymass(idofn) + diagm
     END DO
   CASE (2)
     diagm = diagm/2.d0                    !nodal mass
     DO i=1,3                          !for each element node
       lnode=elem%lnods(i)             !node number
       DO j=1,2
         idofn=iftmp(j,lnode)              !external surface
         IF( idofn > 0 ) ymass(idofn) = ymass(idofn) + diagm
       END DO
     END DO
   CASE (3)
     diagm = diagm*2.d0/3d0                !nodal mass for central node
     val   = diagm/4d0                     !nodal mass for surface nodes
     DO i=1,3                          !for each element node
       lnode=elem%lnods(i)             !node number
       idofn=iftmp(1,lnode)                !mid-surface node
       IF( idofn > 0 ) ymass(idofn) = ymass(idofn) + diagm
       DO j=2,3                            !
         idofn=iftmp(j,lnode)              !external surfaces nodes
         IF( idofn > 0 ) ymass(idofn) = ymass(idofn) + val
       END DO
     END DO
   END SELECT
   elem => elem%next                     !point to next element
 END DO
 RETURN
 END SUBROUTINE tlma24
