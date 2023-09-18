 SUBROUTINE slum24 (nelem,heade,ymass,coord,dtscal)

 ! compute scaled lumped mass vector for element BST-LST (TLF) shell element

 IMPLICIT NONE
 !  dummy arguments
 INTEGER (kind=4), INTENT (IN) :: nelem  !number of elements in the set
 REAL (kind=8), INTENT(IN) :: dtscal
 REAL (kind=8) :: ymass(:,:),      & !nodal mass in each direction
                  coord(:,:)       !nodal coordinates
 TYPE (ele24), POINTER  :: heade     !pointer to first element
 ! local variables
 INTEGER (kind=4) :: ielem,in,jn,kn,isec,osec
 INTEGER (kind=4) :: lnode,idofn,inode
 REAL (kind=8) :: diagm,voli
 REAL (kind=8) :: lmax(4),t1(3),t2(3),young,densi,fac, &
                  a1,a2,b2,poiss,c,r,tcrit,thick,fact  
 TYPE (ele24), POINTER :: elem
 TYPE (section), POINTER :: sec

 fact = 2d0*dtscal**2               !factor
 osec = -1                            !previous material number
 elem => heade                        !point to first element
 DO ielem=1,nelem                     !for each element in the set
   isec = elem%matno                    !associated section
   IF( isec /= osec )THEN               !present and previous material differs
     sec => psecs(isec)%p                 !point to section
     thick = sec%rprop(1)                 !original thickness
     young  = sec%mtbas%prope(1)          !young modulus
     poiss  = sec%mtbas%prope(2)          !poisson ratio
     fac    = young/(1d0-poiss**2)        !plane stress factor
   END IF
   in = elem%lnods(1)                  !first node of the element
   jn = elem%lnods(2)                  !second node of the element
   kn = elem%lnods(3)                  !third node of the element
   t1 = coord(1:3,jn) - coord(1:3,in)  !first side 1-2
   t2 = coord(1:3,kn) - coord(1:3,jn)  !second side 1-3
   a1 = SQRT(DOT_PRODUCT(t1,t1))       !first side length (base length)
   a2 = DOT_PRODUCT(t2,t1)/a1          !proyection of second side over first
   b2 = SQRT(DOT_PRODUCT(t2,t2)-a2**2) !triangle heigth
   ! compute membrane eigenvalues
   c = (b2**2+a2**2+a1*a2+a1**2)/(a1*b2)**2            !center
   r = SQRT(2*a1**3*a2-a1**2*b2**2+a1**4+b2**4        &!radius
              +a2**4+3*poiss**2*a1**2*b2**2           &
              +2*b2**2*a1*a2+3*a2**2*a1**2+2*a2**3*a1 &
              +2*a2**2*b2**2)/(a1*b2)**2
   lmax(1) = (c-r)*fac                           !first eigenvalue
   lmax(2) = (c+r)*fac                           !second eigenvalue
   lmax(3) = c*young/(1d0+poiss)                 !third eigenvalue
   lmax(4) = elem%lb*thick/MIN(a1,b2)            !thickness influence
   tcrit = fact/MAXVAL(lmax(1:3))     !critical time increment
   IF( lmax(4) > 1d0 )tcrit = tcrit/lmax(4)**2
   densi = 1d0/tcrit
   voli = elem%area1*densi*thick         !element mass
   diagm = voli/3.d0                     !nodal mass
   !
   ! ---generates lumped mass matrix proportional to diagonal
   !
   DO inode=1,3                          !for each element node
     lnode=elem%lnods(inode)             !node number
     DO idofn=1,3                        !translational DOFs only
       ymass(idofn,lnode)=ymass(idofn,lnode)+diagm
     END DO
   END DO
   elem => elem%next                     !point to next element
 END DO

 RETURN
 END SUBROUTINE slum24
