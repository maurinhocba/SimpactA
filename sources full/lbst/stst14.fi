 SUBROUTINE stst14( e, iffix, coora, flag, nonrg)
 !-------------------------------------------------------------------
 !
 !     Compute internal variables from stresses for ELEMENT TL CST-BST++
 !
 !-------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: iffix(:)  !boundary conditions
 REAL(kind=8), INTENT(IN)  :: coora(:,:)    !present coordinates
 TYPE (ele14), POINTER  :: e               !pointer to an element
 LOGICAL :: flag !.TRUE.  input stresses are 2nd P-K stresses
                 !.FALSE. input stresses are Kirchhoff stresses
 LOGICAL :: nonrg

 LOGICAL :: sides(3),isidf(3),found
 INTEGER (kind=4) isec,i,j,k,m,n,l,nlayr,ierr,mtype !,nto
 INTEGER (kind=4), SAVE :: osec = -1

 REAL (kind=8) stran(3),stres(3),estra(3),sigma(3),    &
               stra1(6),x(3,6),t(3,2,0:3),t3(3),r1,r2,lb(2),  &
               thnew,dtnew,zk,c1,c2,h(3,3)
 REAL (kind=8), SAVE :: thick,detc,c(4)
 !REAL (kind=8), ALLOCATABLE :: cm(:,:), volfr(:), rr(:,:)
 TYPE(section), POINTER :: sec   !pointer to a section data
 TYPE (mater), POINTER :: mat    !pointer to a material data

 ! Computes layer plastic strains

 isec = e%matno                    !associated material
 IF( isec /= osec )THEN
   osec = isec
   CALL sect_search(isec,found,sec)
   mat => sec%mtbas
   nlayr = sec%iprop(1)              !number of layers
   mtype = mat%mtype                 !type of base material
   IF( mtype == 1)THEN
     c(1:4) = sec%mtbas%prope(7:10)    !elasticity matrix
   ELSE IF( mtype == 5)THEN            !elastic orhthotropic
     ! plane stress elasticity matrix
     c(1) = mat%props(1)
     c(2) = mat%props(2)
     c(3) = mat%props(5)
     c(4) = mat%props(9)
   !ELSE IF( mtype == 9 ) THEN               !mixture theory material
   !  IF( ALLOCATED ( cm ) )DEALLOCATE( cm,volfr,rr )
   !  ALLOCATE( cm(5,nucom),volfr(nucom),rr(5,nucom))
   !  volfr = mat%props(1:nucom)             !volume fractions
   !  CALL vecasi(nucom*5,mat%props(2*nucom+1),cm(1,1)) !elasticity matrix
   !  CALL vecasi(nucom*5,mat%props(7*nucom+1),rr(1,1)) !rotation matrix
   !ELSE IF( mtype == 10 )THEN               ! Elastic Rubber
   !  nto = mat%matdef(6)                    !number of terms in Ogden model
   !  eprop = RESHAPE(mat%prope(7:12),(/2,3/))  ! elastic properties
   !  IF( mat%prope(2) < 0.499 )THEN         ! non- isochoric
   !    isoch = .FALSE.
   !    km = mat%prope(4)*3d0                  ! bulk modulus * 3
   !    DO i=1,nto
   !      aux = eprop(1,i)                  ! exponent alpha
   !      zk  = eprop(2,i)                  ! mu
   !      nr(1,i) = 3d0*zk/km               ! Mu/K
   !      nr(2,i) = aux/3d0                 ! alpha/3
   !      nr(3,i) = 1d0 - aux/3d0           ! 1 - alpha/3
   !      nr(4,i) = 2d0*aux/3d0             !  2*alpha/3
   !      nr(5,i) = 1d0 + 2d0*aux/3d0       ! 1 + 2*alpha/3
   !    END DO
   !  END IF
   END IF
   detc = 1d0/(c(1)*c(3)-c(2)**2)    !to invert constitutive relation
   thick = sec%rprop(1)              !original thickness
 END IF

 c1 = COS(e%angle)   ! COS angle between edge 1-2 and rolling direction
 c2 = SIN(e%angle)   ! SIN angle between edge 1-2 and rolling direction

 !*** computes actual configuration

 isidf = .FALSE.                    !initializes sides boundary conditions
 DO i =1,6                          !for each node in the patch
   n = e%lnods(i)                   !node number
   IF(n /= 0) THEN                  !if node exists
     x(1:3,i) = coora(1:3,n)        !present coordinates
     IF(i > 3)sides(i-3) = .TRUE.   !side element exists
   ELSE
     m = i-3                        !side order
     sides(m) = .FALSE.             !side element does NOT exist
     j = e%lnods(hh(1,m))           !first node of the side
     k = e%lnods(hh(2,m))           !second node of the side
     isidf(m) = iffix(j) == 1 .AND. iffix(k) == 1  !Clamped or S.S. ?
   END IF
 END DO
 CALL stra14(e%a,e%b,x,t3,t,h,sides,isidf,e%ci,nonrg,stra1,e%lb)
 ! Modify curvatures
 stra1(4:6) = stra1(4:6) - e%stra0(1:3)

 thnew = thick*e%lb                  !thickness at new configuration
 dtnew = thnew/nlayr                 !layer thick at new configuration
 zk = -(thnew+dtnew)/2d0             !position for thickness integration

 DO l=1,nlayr                        !for each layer
   zk = zk+dtnew                     !Z coordinate
   stran = stra1(1:3)+stra1(4:6)*zk  !layer U^2
   CALL lgst14(stran,r1,r2,lb,'STST14',ierr)     !Hencky strains, eigen-values/vectors
   stres = e%gausv(1:3,l)            !get stress measures in local frame

   IF( flag ) THEN   !stresses are 2nd Piola-Kirchhoff stresses

     !compute 2nd P-K stresses in the natural Frame
     sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2+2d0*stres(3)*r1*r2
     sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1-2d0*stres(3)*r1*r2
     sigma(3) =(stres(2)-stres(1))*r1*r2+stres(3)*(r1*r1-r2*r2)

     ! Computes Hencky stresses on the natural Frame
     stres(1) = sigma(1)*lb(1)**2
     stres(2) = sigma(2)*lb(2)**2
     IF( ABS(lb(1)-lb(2)) > 1.d-6)THEN
       stres(3) = sigma(3)*(lb(1)**2-lb(2)**2)/2d0/LOG(lb(1)/lb(2))
     ELSE
       stres(3) = sigma(3)*lb(1)*lb(2)
     END IF

     ! Computes Hencky (Kirchhoff) stresses on the Local Frame
     sigma(1) = stres(1)*r1*r1+stres(2)*r2*r2-2d0*stres(3)*r1*r2
     sigma(2) = stres(1)*r2*r2+stres(2)*r1*r1+2d0*stres(3)*r1*r2
     sigma(3) =(stres(1)-stres(2))*r1*r2+stres(3)*(r1*r1-r2*r2)
     stres = sigma
   END IF

   IF( c2 /= 0d0 )THEN   !compute stresses & strains in orthotropic system
     sigma = stres
     stres(1) = sigma(1)*c1*c1+sigma(2)*c2*c2+2d0*sigma(3)*c1*c2
     stres(2) = sigma(1)*c2*c2+sigma(2)*c1*c1-2d0*sigma(3)*c1*c2
     stres(3) =(sigma(2)-sigma(1))*c1*c2+sigma(3)*(c1*c1-c2*c2)
     sigma = stran
     stran(1) = sigma(1)*c1*c1+sigma(2)*c2*c2+2d0*sigma(3)*c1*c2
     stran(2) = sigma(1)*c2*c2+sigma(2)*c1*c1-2d0*sigma(3)*c1*c2
     stran(3) =(sigma(2)-sigma(1))*c1*c2+sigma(3)*(c1*c1-c2*c2)
   END IF
   !compute elastic strains
   estra(1) = (c(3)*stres(1)-c(2)*stres(2))*detc
   estra(2) = (c(1)*stres(2)-c(2)*stres(1))*detc
   estra(3) = stres(3)/c(4)
   e%gausv(1:3,l) = stran - estra        !Plastic strains
 END DO

 RETURN
 END SUBROUTINE stst14
