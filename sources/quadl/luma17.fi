 SUBROUTINE luma17 (heade,emass,sumat,ngaus)

 ! compute lumped mass vector for element 17 (2-D solid TLF)

 IMPLICIT NONE
 !  dummy arguments

 INTEGER (kind=4) :: ngaus
 REAL (kind=8) :: emass(:,:), &    !nodal mass in each direction
                  sumat            !total mass of elements
 TYPE (ele17), POINTER  :: heade   !pointer to first element

 ! local variables

 REAL (kind=8) :: shape(nnode,ngaus),deriv(nnode,2,ngaus), &
                  posgp(2),weigp(2),xita,eta
 INTEGER (kind=4) isec,lnode,inode,i,j,k,osec,ngaud
 REAL (kind=8) :: diagm(nnode),voli,rhoel
 TYPE (ele17), POINTER :: elem
 TYPE (section), POINTER :: sec

 ! -----------------------------
 !ngaud = NINT(SQRT(REAL(ngaus)))
 ngaud = 2
 IF(ngaus == 1)ngaud=1
 CALL gaussq (ngaud ,posgp ,weigp )

 !       gauss points shape and derivatives of nodal functions

 k = 0
 DO i=1,ngaud
   xita = posgp(i)
   DO j=1,ngaud
     k = k+1
     eta = posgp(j)
     CALL shape3(deriv(1,1,k),shape(1,k),xita,eta,nnode)
   END DO
 END DO

 osec = -1
 elem => heade    !point to first element
 DO
   IF (.NOT.ASSOCIATED (elem) ) EXIT     !if all elements processed, exit
   isec = elem%matno                     !associated material
   IF( osec /= isec )THEN
     sec => psecs(isec)%p
     rhoel = sec%mtbas%prope(5)                !material density
     osec = isec
   END IF
   ! ---generates lumped mass matrix proportional to diagonal
   diagm = 0d0
   DO k=1,ngaus
     voli = elem%dvol(k)*rhoel           !Gauss point mass
     diagm = diagm + shape(:,k)*voli
   END DO
   ! ---add to global mass
   DO inode=1,nnode                      !for each element node
     lnode=elem%lnods(inode)             !node number
     emass(1:2,lnode)=emass(1:2,lnode)+diagm(inode)
   END DO
   sumat=sumat+sum(diagm)                !add to total mass
   elem => elem%next                     !point to next element
 END DO
 RETURN
 END SUBROUTINE luma17
