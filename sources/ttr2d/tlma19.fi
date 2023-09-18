 SUBROUTINE tlma19 (heade,ymass,iftmp,ntype,coord)

 ! compute capacity-mass vector for linear triangle element (TLF)

 IMPLICIT NONE
 !  dummy arguments
 INTEGER (kind=4) :: iftmp(:,:), & !solid 2D type problem
                     ntype         !
 REAL (kind=8) :: ymass(:),  &    !lumped nodal mass
                  coord(:,:)      !nodal coordinates (original)
 TYPE (ele19), POINTER  :: heade   !pointer to first element

 ! local variables
 INTEGER (kind=4) isec,i,j,k,osec, &
                  ncurv(3),numct    !number of points in curve
 REAL (kind=8) :: diagm(nnode),voli,rhoel,heatc,val,young,nu,k3,det,r,diag
 TYPE (ele19), POINTER :: elem
 TYPE (section), POINTER :: sec
 TYPE (curve), POINTER :: crvy,crvn,crvh

 osec = -1
 elem => heade    !point to first element
 DO
   IF (.NOT.ASSOCIATED (elem) ) EXIT     !if all elements processed, exit
   isec = elem%matno                    !associated material
   IF( isec /= osec )THEN
     sec => psecs(isec)%p
     young = sec%mtbas%prope(1)          !young modulus
     nu    = sec%mtbas%prope(2)          !nu
     rhoel = sec%mtbas%prope(5)          !material density
     heatc = sec%mtbas%prope(7)          !material capacity
     k3    = 3d0*sec%mtbas%prope(4)      !bulk modulus * 3
     !associated curves for thermal dependant properties
     ncurv(1) = INT(sec%mtbas%props(1))     !young modulus
     ncurv(2) = INT(sec%mtbas%props(2))     !poisson ratio
     ncurv(3) = INT(sec%mtbas%props(7))     !heat capacity
     CALL cur_point(ncurv(1),sec%mtbas%chead,crvy) !pointer to Young curve
     CALL cur_point(ncurv(2),sec%mtbas%chead,crvn) !pointer to Poisson curve
     CALL cur_point(ncurv(3),sec%mtbas%chead,crvh) !pointer to Heat Capac. curve
     osec = isec
   END IF
   numct = SUM(ncurv)
   IF( numct > 0 )THEN !compute properties as a function of temperature
     ! temperature at element center
     val = SUM(tempe(1,elem%lnods(1:3)))/3d0
     i = 1
     IF(ncurv(1) > 0) young = inte_cr(crvy%val,crvy%np,val,i)
     IF(ncurv(2) > 0) nu    = inte_cr(crvn%val,crvn%np,val,i)
     IF(ncurv(3) > 0) heatc = inte_cr(crvh%val,crvh%np,val,i)
     IF(ANY(ncurv(1:2) > 0)) k3 = 3d0*(young/(1d0+nu)/2d0) !bulk modulus * 3
   END IF
   !voli = elem%area1*rhoel/det    !total element mass (from element 20)
   det = EXP((elem%stint(1)+elem%stint(2)+elem%stint(4))/k3)
   voli=(elem%area1*rhoel/det)*heatc     !element capacity
   !
   ! ---generates lumped mass matrix proportional to diagonal
   IF( ntype == 3 ) THEN  !for axilsymmetric problems
     diagm = 0d0
     diag = voli/SUM(coord(1,elem%lnods(1:3)))/4d0
     DO i=1,3     !for each side
       j = elem%lnods(kk(1,i))
       k = elem%lnods(kk(2,i))
       r = coord(1,j) + coord(1,k)
       diagm(kk(1,i)) = diagm(kk(1,i)) + diag*r
       diagm(kk(2,i)) = diagm(kk(2,i)) + diag*r
     END DO
   ELSE
     diagm = voli/3d0
   END IF
   !
   ! ---add to global mass
   DO i=1,nnode                      !for each element node
     j=elem%lnods(i)             !node number
     k=iftmp(1,j)
     IF( k > 0 ) ymass(k)=ymass(k)+diagm(i)
   END DO
   elem => elem%next                     !point to next element
 END DO
 RETURN
 END SUBROUTINE tlma19
