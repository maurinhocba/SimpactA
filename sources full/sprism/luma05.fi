 SUBROUTINE luma05 (heade,ymass,sumat)

 ! compute lumped mass vector for element 05 (3-D solid TLF)
 ! use standard 2x1   integration rule

 IMPLICIT NONE
 !  dummy arguments

 REAL (kind=8) :: ymass(:,:), &    !nodal mass in each direction
                  sumat            !total mass of elements
 TYPE (ele05), POINTER  :: heade   !pointer to first element
 ! local variables
 INTEGER (kind=4), PARAMETER ::  ngaus=2
 REAL (kind=8) :: shape(nnb,ngaus),deriv(nnb,3), &
                  xita,eta,posgp(ngaus),weigp(ngaus),sh(3,ngaus)
 INTEGER (kind=4) isec,lnode,inode,i,k,osec
 REAL (kind=8) :: diagm(nnb),voli,dens
 TYPE (ele05), POINTER :: elem
 TYPE(section), POINTER :: sec
 LOGICAL :: newm

 ! -----------------------------

 CALL gaussq (ngaus ,posgp ,weigp )

 !       gauss points shape and derivatives of nodal functions

 xita = 1d0/3d0
 eta  = 1d0/3d0
 weigp = weigp/2d0
 DO i=1,ngaus
    sh(1,i) = -0.5d0*(1d0-posgp(i))*posgp(i)
    sh(2,i) =  1d0-posgp(i)**2
    sh(3,i) =  0.5d0*(1d0+posgp(i))*posgp(i)
    CALL shape4(nnb,shape(1,i),deriv(1,1),xita,eta,posgp(i))
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
   ! ---generates lumped mass matrix proportional to diagonal
   diagm = 0d0
   DO k=1,ngaus
     voli = DOT_PRODUCT(elem%dvol,sh(:,k))*weigp(k)*dens   !GP mass
     diagm = diagm + shape(:,k)*voli
   END DO
   ! ---add to global mass
   DO inode=1,nnb                          !for each element node
     lnode = elem%lnods(inode)             !node number
     ymass(1:3,lnode) = ymass(1:3,lnode) + diagm(inode)
   END DO
   sumat = sumat + SUM(diagm)            !add to total mass
   elem => elem%next                     !point to next element
 END DO
 RETURN
 END SUBROUTINE luma05
