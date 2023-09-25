 SUBROUTINE masm18 (heade,emass,mass,sumat,coord,lumpd)

 ! compute mass vector/matrix for element 18 (3-D solid TLF)
 ! use standard 2x2x2 integration rule

 IMPLICIT NONE
 !  dummy arguments
 LOGICAL, INTENT(IN) :: lumpd
 REAL (kind=8) :: emass(:,:), &    !nodal mass in each direction
                  mass(:),    &    !assembled consistent mass matrix
                  coord(:,:), &    !nodal coordinates
                  sumat            !total mass of elements
 TYPE (ele18), POINTER  :: heade   !pointer to first element
 ! local variables
 INTEGER (kind=4), PARAMETER :: ngaud=2, ngaus=8, nnode=8, ndime=3
 REAL (kind=8) :: shape(nnode,ngaus),deriv(nnode,3,ngaus),wei(ngaus), &
                  posgp(ngaud),weigp(ngaud),xita,eta,zeta,cd(nnode,3),x(3,nnode)
 INTEGER (kind=4) isec,lnode,inode,i,ii,j,k,l,n,osec,err
 REAL (kind=8) :: diagm(nnode),voli,dens
 REAL (kind=8), ALLOCATABLE :: a(:),ms(:)
 TYPE (ele18), POINTER :: elem
 TYPE(section), POINTER :: sec
 LOGICAL :: newm
 INTEGER (kind= 4) posmt
 posmt(k,i,j) = (k-i)*(i-1)/2+j  !position i,j in mass(k/2 x k/2)

 ! -----------------------------

 CALL gaussq (ngaud ,posgp ,weigp )

 !       gauss points shape and derivatives of nodal functions

 k = 0
 DO l=1,ngaud
   zeta = posgp(l)
   DO j=1,ngaud
     eta = posgp(j)
     DO i=1,ngaud
       xita = posgp(i)
       k = k+1
       CALL shape5(deriv(1,1,k),shape(1,k),xita,eta,zeta,nnode)
       wei(k) = weigp(i)*weigp(j)*weigp(l)
     END DO
   END DO
 END DO

 n = nnode*(nnode+1)/2                    !size of array a
 l = (nnode*ndime*(nnode*ndime+1))/2      !size of mass matrix
 ALLOCATE(a(n),ms(l))

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
     ii = 0                   !Initializes Position in matrix
     DO n = 1,nnode           !For each node
       IF(.NOT.lumpd) THEN    !Consistent matrix
           DO l = n,nnode       !for each node
             ii = ii+1          !position in matrix
             a(ii) = a(ii) + shape(n,k)*shape(l,k)*voli
           END DO
        END IF     !.NOT.lumpd
     END DO
   END DO
   ! ---add to global mass
   DO inode=1,nnode                        !for each element node
     lnode = elem%lnods(inode)             !node number
     emass(1:3,lnode) = emass(1:3,lnode) + diagm(inode)
   END DO

   IF(.NOT.lumpd) THEN
                 !consistent mass matrix  Assembles over MASS(maxa)
     i = 0                            ! Position in matrix A initialization
     k = 0                            ! elemental DOF initialization
     DO n = 1,nnode                   ! For each nodee
       j = (n-1)*3                ! Initializes Column position in M
       DO l = n,nnode                 ! For each node
         i = i+1                      ! Matrix position in A
         k = (l-1)*3              ! Initializes Row position in M
         DO ii = 1,3              ! For each direction in space
           ms(posmt(48,j+ii,k+ii)) = a(i)
         END DO
       END DO
     END DO
     CALL masaux(8,elem%lnods,24,3,ms(1),mass(1) )
   END IF

   sumat = sumat + sum(diagm)            !add to total mass
   elem => elem%next                     !point to next element
 END DO
 RETURN
 END SUBROUTINE masm18
