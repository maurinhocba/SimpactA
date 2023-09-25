 SUBROUTINE luma19 (heade,ymass,sumat,ntype,coord)

 ! compute lumped mass vector for element 19 (2-D solid TLF)

 IMPLICIT NONE
 !  dummy arguments

 INTEGER (kind=4), INTENT(IN) :: ntype
 REAL (kind=8) :: ymass(:,:), &  !nodal mass in each direction
                  sumat,      &  !total mass of elements
                  coord(:,:)     !nodal coordinates (original)
 TYPE (ele19), POINTER  :: heade   !pointer to first element

 ! local variables
 INTEGER (kind=4) isec,i,j,k,osec
 REAL (kind=8) :: diagm(nnode),voli,rhoel,r,diag,setmas,det,k3
 TYPE (ele19), POINTER :: elem
 TYPE (section), POINTER :: sec

 ! -----------------------------

 osec = -1
 setmas = sumat
 elem => heade    !point to first element
 DO
   IF (.NOT.ASSOCIATED (elem) ) EXIT     !if all elements processed, exit
   isec = elem%matno                     !associated material
   IF( osec /= isec )THEN
     sec => psecs(isec)%p
     rhoel = sec%mtbas%prope(5)                !material density
     k3    = 3d0*sec%mtbas%prope(4)            !bulk modulus
     osec = isec
   END IF
   det = EXP((elem%stint(1)+elem%stint(2)+elem%stint(4))/k3)
   voli = elem%area1*rhoel/det             !total mass
   sumat = sumat + voli                    !set mass
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
   ! ---add to global mass
   DO i=1,nnode                      !for each element node
     j=elem%lnods(i)             !node number
     ymass(1:2,j)=ymass(1:2,j)+diagm(i)
   END DO
   elem => elem%next                     !point to next element
 END DO
 setmas = sumat - setmas
 WRITE(55,"(' total mass at this moment: ',e15.6)")setmas
 RETURN
 END SUBROUTINE luma19
