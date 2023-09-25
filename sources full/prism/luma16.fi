 SUBROUTINE luma16 (heade,ymass,sumat,coord)

 ! compute lumped mass vector for element 16 (3-D solid TLF)
 ! use standard 2x2x2 integration rule

 IMPLICIT NONE
 !  dummy arguments

 REAL (kind=8) :: ymass(:,:), &    !nodal mass in each direction
                  coord(:,:), &    !nodal coordinates
                  sumat            !total mass of elements
 TYPE (ele16), POINTER  :: heade   !pointer to first element
 ! local variables
 INTEGER (kind=4), PARAMETER :: ngaud = 2, ngaus=2
 REAL (kind=8) :: shape(nnode,ngaus),deriv(nnode,3,ngaus),wei(ngaus), &
                  posgp(ngaud),weigp(ngaud),xita,eta,zeta,cd(nnode,3),x(3,nnode)
 INTEGER (kind=4) isec,lnode,inode,i,k,osec,err
 REAL (kind=8) :: diagm(nnode),voli,dens
 TYPE (ele16), POINTER :: elem
 TYPE(section), POINTER :: sec
 LOGICAL :: newm

 ! -----------------------------

 CALL gaussq (ngaud ,posgp ,weigp )

 !       gauss points shape and derivatives of nodal functions

 xita = 1d0/3d0
 eta  = 1d0/3d0
 DO i=1,ngaud
    zeta = posgp(i)
    CALL shape4(nnode,shape(1,i),deriv(1,1,i),xita,eta,zeta)
    wei(i) = 0.5d0
 END DO

 osec = -1
 elem => heade    !point to first element
 DO
   IF (.NOT.ASSOCIATED (elem) ) EXIT     !if all elements processed, exit
   isec = elem%matno                     !associated material
   newm = isec /= osec
   IF ( newm )THEN
     sec => psecs(isec)%p
     IF( sec%mabas == -1 )THEN
       dens = sec%rprop(3)
     ELSE
       dens = sec%mtbas%prope(5)
     END IF
     osec = isec
   END IF
   !
   x = coord(:,elem%lnods)
   ! ---generates lumped mass matrix proportional to diagonal
   diagm = 0d0
   DO k=1,ngaus
     CALL jacob5(cd(1,1),deriv(1,1,k),voli,x(1,1),nnode,999,err)
     voli = voli*dens*wei(k)           !Gauss point mass
     diagm = diagm + shape(:,k)*voli
   END DO
   ! ---add to global mass
   DO inode=1,nnode                        !for each element node
     lnode = elem%lnods(inode)             !node number
     ymass(1:3,lnode) = ymass(1:3,lnode) + diagm(inode)
   END DO
   sumat = sumat + sum(diagm)            !add to total mass
   elem => elem%next                     !point to next element
 END DO
 RETURN
 END SUBROUTINE luma16
