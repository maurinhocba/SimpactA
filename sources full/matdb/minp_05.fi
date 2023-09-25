 SUBROUTINE minp_05 ( mat)
 !            --------------------------------
 !            input material properties for
 !            elastic - plastic, orthotropic
 !            elastic stress-strain curve: Linear
 !            plastic stress-strain curve: linear, exponential, saturation
 !            --------------------------------
 ! MATDEF: 1: 2,3    2-D, 3-D
 !         2: 3      Orthotropic
 !         3: 1,2,3  Elastic, von Mises, Hill-48
 !         4: 1,2,3  No Hard, Linear, Ludwik-Nadai
 !         5: 1      No kinematic Hardening
 !         6: 1,2    Elastic, elastic-plastic
 !         7: 0
 !         8: 1234   submodel
 ! PROPE: 1:Max Young or 0D0       5:Density   6:Alpha
 !        7-12 : E1, E2, E3, NU12, NU13, NU23 or  D1111, D1122, D2222, D1212, D3333, D1133,
 !        13-15: G12, G13, G23                or  D2233, D1313, D2323,
 !        Elasticity matrices,  lower part only has non zero values
 !       16-19: 2-D plane stress C11 C12 C22 C33
 !       20-28: 2-D plane stress (3x3)
 !       29-44: 2-D Plane strain or axi-symmetric (4x4)
 !       45-80: 3-D (6x6)
 ! PROPP:  1-4: isotropic hardening parameters, according to MATDEF(4)
 !             (1) UNIAX
 !             (2) UNIAX, HARDS
 !             (3) Cons0, Efref, Consn
 !             (4) UNIAX, HARDS, SATEXP, REFLST
 !           5: death time
 !         6-7: kinematic hardening parameters
 !        8-11: s90,biuni,sxy,s45
 !       12-15: rm,r0,r45,r90
 !       17:21: Isotropic hardening parameters
 !       22-24: Strength ratios
 !       25-29: AE matrix
 ! PROPS: Nothing
 IMPLICIT NONE
 ! dummy arguments
 TYPE (mater), POINTER :: mat
 ! local variables
 REAL (kind=8) :: cm1,cm2,cm3,aux

 ! check Material Definitions

 IF(mat%matdef(2) == 0 ) mat%matdef(2) = 3   !orthotropic
 IF(mat%matdef(2) /= 3 )CALL runend('MINP05: mtype = 5 must be ORTHOTROPIC ')
 IF(mat%matdef(3)  > 3 )CALL runend('MINP05: mtype = 5 cannot use HILL79/90')
 IF(mat%matdef(4)  > 4 )CALL runend('MINP05: mtype = 5 cannot use IHPOINT  ')
 IF(mat%matdef(5)  > 1 )CALL runend('MINP05: mtype = 5 cannot use KIN-HARD ')

 ! get memory for constants & read elastic properties
 mat%matdef(9)  = 80                        !size of array PROPE
 ALLOCATE( mat%prope(mat%matdef(9)) )       !get memory for PROPE
 mat%prope = 0d0                            !initializes PROPE
 CALL mat_ela(mat%prope,mat%matdef)         !read elastic constants
 CALL setm_05(mat%prope(7:15),mat%prope(20:80),mat%matdef(8)) !compute elastic matrices

 ! repeat plane stress matrix for space mapping plasticity
 IF( mat%matdef(8) == 1 )THEN         !plane stress elasticity matrix was READ
   mat%prope(16:19) =    mat%prope(7:10)
 ELSE
   mat%prope(16) = mat%prope(20)      !C(1,1)
   mat%prope(17) = mat%prope(21)      !C(2,1)
   mat%prope(18) = mat%prope(24)      !C(2,2)
   mat%prope(19) = mat%prope(28)      !C(3,3)
 END IF
 ! get memory for constants & read plastic properties

 IF( mat%matdef(3) /= 1 )THEN !plastic
   mat%matdef(10)  = 29                     !size of array PROPP
   ALLOCATE( mat%propp(mat%matdef(10)))     !get memory for PROPP
   mat%propp = 0d0                          !initializes PROPP
   CALL mat_pla(mat%propp,mat%matdef)       !read plastic constants
   ! next values to use space mapping plasticity
   mat%propp(17:20) = mat%propp(1:4)        !repeat isotropic hardening data
   mat%propp(21) = REAL(mat%matdef(4))      !isotropic hardening model
   mat%propp(22) = 1d0                      !S0
   mat%propp(23) = mat%propp(8)             !S90
   mat%propp(24) = mat%propp(10)            !S45
   ! Isotropic C
   mat%prope(20) = (mat%prope(16)-mat%prope(17))/2d0   !shear term
   ! and its inverse
   aux = mat%prope(16)*mat%prope(16)-mat%prope(17)*mat%prope(17) !det
   cm1 =  mat%prope(16)/aux
   cm2 = -mat%prope(17)/aux
   cm3 =  1d0/mat%prope(20)
   ! AE
   mat%propp(25) =  cm1*mat%prope(16)/mat%propp(22) + cm2*mat%prope(17)/mat%propp(23)
   mat%propp(26) =  cm1*mat%prope(17)/mat%propp(22) + cm2*mat%prope(18)/mat%propp(23)
   mat%propp(27) =  cm2*mat%prope(16)/mat%propp(22) + cm1*mat%prope(17)/mat%propp(23)
   mat%propp(28) =  cm2*mat%prope(17)/mat%propp(22) + cm1*mat%prope(18)/mat%propp(23)
   mat%propp(29) =  mat%prope(19)*cm3/mat%propp(22)
 ELSE                         !elastic
   mat%matdef(10)  = 1                       !set size to 1
   ALLOCATE( mat%propp(mat%matdef(10)) )
 END IF
 !  reserve space for (not used) array props
 mat%matdef(11) = 1                          !size of array PROPS
 ALLOCATE( mat%props(mat%matdef(11)) )       !get memory for PROPS

 RETURN
 END SUBROUTINE minp_05

 SUBROUTINE setm_05(epr,prop,model)
 !            --------------------------------
 !            input material properties for
 !            plane stress/strain, axisymmetric
 !            linear elastic orthotropic
 !            -- material can be used for:
 !               a) 2d continuum elements
 !               b) 3d shell elements
 !               c) composite laminates
 !
 IMPLICIT NONE
 REAL (kind=8), INTENT(IN OUT) :: epr(:)
 REAL (kind=8), INTENT(OUT) :: prop(:)
 INTEGER (kind=4), INTENT(IN) :: model

 LOGICAL allparams

 ! models
 ! 1 - plane stress
 ! 2 - plane strain
 ! 3 - axisymmetric
 ! 4 - 3D

 ! 3-D form

 IF (model == 0 ) THEN

   ! 3D const matrix if possible (if all the parameters input)

   ! if all parameters defined (3D)
   allparams = (epr(1)*epr(2)*epr(3)*epr(7) > 1.e-8)   ! 3D computed if exist E1, E2 and E3 - XMG

   IF( allparams ) THEN
     CALL omat3d(prop(26:61),epr(1:9),allparams,.TRUE.)       !computes 3-D matrix
     CALL omat2d2(prop(26:61),prop(1:9),prop(10:25))
   ELSE
    ! 2D forms
     CALL omat2d1(epr(1),epr(2),epr(4),epr(7),prop(1:9))
     CALL omat3d(prop(26:61),epr(1:9),allparams,.TRUE.)         !computes 3-D matrix
   END IF

 ELSE
   ! 2D forms
   CALL omat2d3(epr(1:9),prop(1:9),prop(10:25))  !transfer user defined D matrix

 END IF

 RETURN
 END SUBROUTINE setm_05

 SUBROUTINE omat3d(mtx,epr,all,prt)
 !************************************************************************
 !
 !****this routine sets up the elastic matrix for 3-d problems
 !
 !************************************************************************
 IMPLICIT NONE
 REAL (kind=8), INTENT(OUT) :: mtx(36)
 REAL (kind=8), INTENT(IN) :: epr(9)
 LOGICAL, INTENT(IN) :: all,prt

 REAL (kind=8) exx, pxy, eyy, pxz, ezz, pyz, gxy, gxz, gyz, pyx,  pzx, pzy, &
               delta
 LOGICAL :: std

 LOGICAL :: problems(7)
 !REAL (kind=8) d(3,3)
 !***initialize the elastic matrix

 mtx = 0d0         !initializes all values

 !   convention used      p_ij = -eps_j/eps_i
 !                        p_ji*e_ii = p_ij*e_jj

 ! get data
 exx = epr(1)   ! E1
 pxy = epr(4)   ! v12
 eyy = epr(2)   ! E2
 pxz = epr(5)   ! v12
 ezz = epr(3)   ! E3
 pyz = epr(6)   ! v23
 gxy = epr(7)   ! G12
 gxz = epr(8)   ! G13
 gyz = epr(9)   ! G23
 IF( .NOT.all )THEN
   IF( eyy == 0d0 )RETURN  !unidirectional model
   std = eyy < exx     !x is the fiber direction
   IF( std )THEN
     pxz = pxy      ! v13
     ezz = eyy      ! E3
     pyz = 0.25     ! v23
     gxz = gxy      ! G13
     gyz = 0.4*eyy  ! G23
   ELSE
     pyz = pxy*(eyy/exx) ! v23
     ezz = exx           ! E3
     pxz = 0.25          ! v13
     gyz = gxy           ! G23
     gxz = 0.4*exx       ! G13
   END IF
 END IF
 ! compute others poisson's ratios
 pyx = pxy*(eyy/exx)   !v21
 pzx = pxz*(ezz/exx)   !v31
 pzy = pyz*(ezz/eyy)   !v32
 ! checks
 delta = 1d0 - pxy*pyx-pyz*pzy-pxz*pzx-pyx*pxz*pzy-pxy*pyz*pzx
 problems(1) = exx <= eyy*pxy**2
 problems(2) = eyy <= ezz*pyz**2
 problems(3) = ezz <= exx*pzx**2
 problems(4) = delta <= 0d0
 problems(5) = 1d0 <= pxy*pyx
 problems(6) = 1d0 <= pzy*pyz
 problems(7) = 1d0 <= pxz*pzx
 IF( ANY(problems) ) THEN
   WRITE( lures, "('WARNING: Material may not satisfy stability! ')")
   IF(problems(1)) WRITE( lures, "(' E_11 <= E_22 * NU_12^2  , NU_12=',f7.3)")pxy
   IF(problems(2)) WRITE( lures, "(' E_22 <= E_33 * NU_23^2  , NU_23=',f7.3)")pyz
   IF(problems(3)) WRITE( lures, "(' E_33 <= E_11 * NU_31^2  , NU_31=',f7.3)")pzx
   IF(problems(4)) WRITE( lures, "(' Delta <= 0d0  , delta=',f7.3)")delta
   IF(problems(5)) WRITE( lures, "(' NU_12 * NU_21 >= 1, NU_12=',f7.3,' NU_21=',f7.3)")pxy,pyx
   IF(problems(6)) WRITE( lures, "(' NU_23 * NU_32 >= 1, NU_23=',f7.3,' NU_32=',f7.3)")pyx,pzy
   IF(problems(7)) WRITE( lures, "(' NU_31 * NU_13 >= 1, NU_31=',f7.3,' NU_13=',f7.3)")pzx,pxz
 END IF
! !this is a possible way to compute elasticity matrix
!  d(1:3,1) = (/  1d0/exx,-pxy/exx,-pxz/exx/)
!  d(1:3,2) = (/ -pyx/eyy, 1d0/eyy,-pyz/eyy/)
!  d(1:3,3) = (/ -pzx/ezz,-pzy/ezz, 1d0/ezz/)
!  CALL posin1 (d,3,3)
! !or explicitly
! order 11-22-33-12-13-23  (this is the way it is used)
 mtx(1) = exx*(1d0 - pyz*pzy)/delta  !(1,1)   11-11
 mtx(2) = eyy*(pxy + pxz*pzy)/delta  !(2,1)   11-22
 mtx(3) = ezz*(pxz + pxy*pyz)/delta  !(3,1)   11-33
 mtx(7) = mtx(2)                     !(1,2)   22-11
 mtx(8) = eyy*(1d0 - pxz*pzx)/delta  !(2,2)   22-22
 mtx(9) = ezz*(pyz + pxz*pyx)/delta  !(3,2)   22-33
 mtx(13)= mtx(3)                     !(1,3)   33-11
 mtx(14)= mtx(9)                     !(2,3)   33-22
 mtx(15)= ezz*(1d0 - pxy*pyx)/delta  !(3,3)   33-33
 mtx(22)= gxy                        !(4,4)   12-12
 mtx(29)= gxz                        !(5,5)   13-13
 mtx(36)= gyz                        !(6,6)   23-23
 IF( prt )THEN
   WRITE(lures,"(' 3D Elasticity Matrix ')")
   WRITE(lures,"(3e15.6,/,3e15.6,/,3e15.6)")mtx(1:3),mtx(7:9),mtx(13:15)
   WRITE(lures,"(45x,e15.6,/,60x,e15.6,/,75x,e15.6)")mtx(22),mtx(29),mtx(36)
 END IF
 RETURN
 END SUBROUTINE omat3d

 SUBROUTINE omat2d1 (e1,e2,nu12,g12,c)
 IMPLICIT NONE
 !
 !       stress-strain law for plane stress orthotropic materials
 !       in material axes
 !
 REAL (kind=8) :: c(9),e1,e2,nu12,g12

 REAL (kind=8) :: nu21,delta,aux
 !
 c = 0

 nu21=nu12*e2/e1
 IF( nu21 > 0.5D0)THEN
   WRITE(*,"('WARNING NU21 > 0.5,   NU_21= ',f6.3)")NU21
   WRITE(lures,"('WARNING NU21 > 0.5,   NU_21= ',f6.3)")NU21
 END IF
 aux=nu21*e1
 delta = 1.d0 - nu12*nu21
 IF( delta < 0.05D0)THEN
   WRITE(*,"('WARNING 1.d0 - nu12*nu21 < 0.05   DELTA=',f10.4)")delta
   WRITE(lures,"('WARNING 1.d0 - nu12*nu21 < 0.05   DELTA=',f10.4)")delta
   IF( delta <= 0d0 ) STOP
 END IF

 c(1) = e1/delta      ! c(1,1)
 c(2) = aux/delta     ! c(2,1)
 c(4) = c(2)          ! c(1,2)
 c(5) = e2/delta      ! c(2,2)
 c(9) = g12           ! c(3,3)

 WRITE(lures,"(' 2D Plane-Stress Elasticity Matrix ')")
 WRITE(lures,"(2e15.6,/,2e15.6)")c(1:2),c(4:5)
 WRITE(lures,"(30x,e15.6)")c(9)
 RETURN
 END SUBROUTINE omat2d1

 SUBROUTINE omat2d2 (c6,c3,c4)
 IMPLICIT NONE
 !
 !       stress-strain law for plane strain orthotropic materials
 !       in material axes
 !
 REAL (kind=8) c6(36),c4(16),c3(9)
 REAL (kind=8) a

 c4 = 0d0

 !       directly extract Plain strain matrix from 3D
 c4( 1) = c6( 1)                 !(1,1)
 c4( 6) = c6( 8)                 !(2,2)
 c4(11) = c6(22)                 !(3,3)
 c4(16) = c6(15)                 !(4,4)
 c4( 2) = c6( 2)                 !(2,1)
 c4( 4) = c6( 3)                 !(4,1)
 c4( 8) = c6( 9)                 !(4,2)
 c4( 5) = c4( 2)                 !(1,2)
 c4(13) = c4( 4)                 !(1,4)
 c4(14) = c4( 8)                 !(2,4)
 WRITE(lures,"(' 2D Plane-Strain Elasticity Matrix ')")
 WRITE(lures,"(4e15.6)")c4(1:16)
 !
 !     for plane stress analysis condense stress-strain matrix (S_33 = 0)
 !
 a = c4(4)/c4(16)                  !c4(1,4)/c4(4,4)
 c3(1:3) = c4(1:3) - c4(13:15)*a   !c4(1,1:3) - c4(4,1:3)*a
 a = c4(8)/c4(16)                  !c4(2,4)/c4(4,4)
 c3(4:6) = c4(5:7) - c4(13:15)*a   !c4(2,1:3) - c4(4,1:3)*a
 a = c4(12)/c4(16)                 !c4(3,4)/c4(4,4)
 c3(7:9) = c4(9:11) - c4(13:15)*a  !c4(3,1:3) - c4(4,1:3)*a
 WRITE(lures,"(' 2D Plane-Stress Elasticity Matrix ')")
 WRITE(lures,"(3e15.6)")c3(1:9)

 RETURN
 !
 END SUBROUTINE omat2d2

 SUBROUTINE omat2d3 (dijkl,c3,c4)
 IMPLICIT NONE
 !
 !       stress-strain law for plane strain orthotropic materials
 !       in material axes
 !
 REAL (kind=8) :: c3(:),c4(:),dijkl(:)

 !     dijkl(1)   d1111
 !     dijkl(4)   d1122
 !     dijkl(2)   d2222
 !     dijkl(5)   d1133
 !     dijkl(6)   d2233
 !     dijkl(3)   d3333
 !     dijkl(7)   d1212

 REAL (kind=8) :: a

 c3 = 0d0  !initializes both matrices
 c4 = 0d0

 c4( 1) = dijkl(1)  !(1,1)    !transfer user defined matrix
 c4( 2) = dijkl(4)  !(2,1)
 c4( 4) = dijkl(5)  !(4,1)
 c4( 5) = c4(2)     !(2,1)
 c4( 6) = dijkl(2)  !(2,2)
 c4( 8) = dijkl(6)  !(4,2)
 c4(11) = dijkl(7)  !(3,3)
 c4(13) = c4(4)     !(1,4)
 c4(14) = c4(8)     !(4,2)
 c4(16) = dijkl(3)  !(4,4)
 WRITE(lures,"(' 2D Plane-Strain Elasticity Matrix ')")
 WRITE(lures,"(4e15.6)")c4(1:16)
 !
 !     for plane stress analysis condense stress-strain matrix
 !
 a = c4(4)/c4(16)                  !c4(1,4)/c4(4,4)
 c3(1:3) = c4(1:3) - c4(13:15)*a   !c4(1,1:3) - c4(4,1:3)*a
 a = c4(8)/c4(16)                  !c4(2,4)/c4(4,4)
 c3(4:6) = c4(5:7) - c4(13:15)*a   !c4(2,1:3) - c4(4,1:3)*a
 a = c4(12)/c4(16)                 !c4(3,4)/c4(4,4)
 c3(7:9) = c4(9:11) - c4(13:15)*a  !c4(3,1:3) - c4(4,1:3)*a
 WRITE(lures,"(' 2D Plane-Stress Elasticity Matrix ')")
 WRITE(lures,"(3e15.6)")c3(1:9)

 RETURN
 END SUBROUTINE omat2d3
