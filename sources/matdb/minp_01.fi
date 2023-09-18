 SUBROUTINE minp_01 ( mat )
 !            --------------------------------
 !            input material properties for
 !            elastic isotropic
 !            plastic orthotropic:  von Mises - Hill 48
 !            elastic stress-strain curve: linear
 !            isotropic hardening: Linear - Ludwik Nadai, Linear+Saturation
 !            no-kinematic hardening yet
 !            --------------------------------
 ! MATDEF: 1: 3      3-D
 !         2: 1,2    Rigid,Isotropic
 !         3: 1--6   Elastic, von Mises, Hill-48, Hill-79, Hill-90, BBC'03
 !         4: 1,2,3,4  No Hard, Linear, Ludwik-Nadai, linear + saturation
 !            5,6,7    Points, Holomon , Voce  (this 3 Not valid for Shear Def. shells)
 !         5: 1,2      Linear kinematic Hardening    (valid only for Shear Def shells)
 !         6: 1,2,3    Elastic, Elastic-Plastic, visco-plastic
 !         7: 0
 !         8: 0,1    Associative, Non-associative
 ! PROPE: 1:Young 2:Poiss 3:Shear G  4:Bulk K  5:Density   6:Alpha
 !        7-10: plane stress C (C11,C12,C22,C33)
 ! PROPP:  1-4: isotropic hardening parameters, according to MATDEF(4)
 !             (1)IHNONE: Uniax
 !             (2)IHLINE: Uniax, Hards
 !             (3)IHLUDW: Cons0, Efref, Consn
 !             (4)IHSATU: Uniax, Hards, SatExp, ReFlSt
 !             (5)IHPOIN: Uniax, data is included in CHEAD
 !             (6)IHHOLO: Cons0,     0, ConsN
 !             (7)IHVOCE: Cons0, Uniaz, SatEx
 !           5: death time
 !         6-7: kinematic hardening parameters
 !        8-11: s90,biuni,sxy,s45: strength ratios
 !       12-15: rm,r0,r45,r90    : Lankford parameters
 !       16-21: Hill coefficients for yield function
 !       22-27: Hill coefficients for potential function
 IMPLICIT NONE
 ! dummy arguments
 TYPE (mater), POINTER :: mat
 ! local variables
 REAL(kind=8) :: nu2

 ! check Material Definitions

 IF(mat%matdef(2) == 0 ) mat%matdef(2) = 2  !isotropic
 IF(mat%matdef(3) == 0 ) mat%matdef(3) = 1  !elastic
 IF( mat%matdef(2) == 1 )THEN               !for a rigid material
   mat%matdef(9)  = 5                     !size of array PROPE (only the fifth value is needed)
 ELSE                                       !for a deformable material
   IF(mat%matdef(2) /= 2 )CALL runend('MINP01: mtype = 1 must be ISOTROPIC')
   IF(mat%matdef(3)  > 6 )CALL runend('MINP01: mtype = 1 invalid Yield Function') !FR
   IF(mat%matdef(4)  > 7 )CALL runend('MINP01: mtype = 1 invalid Iso-Hard law')
   IF(mat%matdef(5)  > 3 )CALL runend('MINP01: mtype = 1 invalid Kin-Hard law')
   mat%matdef(9)  = 10                    !size of array PROPE
 END IF
 ! get memory for constants & read elastic properties
 ALLOCATE( mat%prope(mat%matdef(9)) )   !get memory for PROPE
 mat%prope = 0d0                        !initializes PROPE
 CALL mat_ela(mat%prope,mat%matdef)     !read elastic constants
 IF( mat%matdef(2) /= 1 )THEN ! not rigid material
   ! compute plane stress elasticity matrix
   nu2 = 1d0-mat%prope(2)*mat%prope(2)                   !1-nu^2
   mat%prope(7)  = mat%prope(1)/nu2                      !E/(1-nu^2)
   mat%prope(8)  = mat%prope(7)*mat%prope(2)             !*nu
   mat%prope(9)  = mat%prope(7)                          !==C11
   mat%prope(10) = mat%prope(1)/(2d0*(1d0+mat%prope(2))) !E/2(1+nu)
   WRITE(58,"('label :',i5,/,5x,'E =',e15.5,/,5x,'nu=',f10.5,/,5x,'G =',e15.5,/,5x,'K =',e15.5)")mat%matno,mat%prope(1:4)
   nu2 = (mat%prope(7)*mat%prope(9)-mat%prope(8)**2)*mat%prope(10)         !det
   WRITE(58,"('   Plane Stress Matrix',/,5x,2e15.5,/,20x,E15.5,/,35x,E15.5,/,'   Det:',e15.5)")mat%prope(7:10),nu2
 END IF
   ! get memory for constants & read plastic properties
 IF( mat%matdef(3) > 1 )THEN    !if not ELASTIC
   mat%matdef(10)  = 27                  !size of array PROPP
   IF(mat%matdef(6) == 3) mat%matdef(10) = 30 !add space for visco-plast prop
   ALLOCATE( mat%propp(mat%matdef(10)) ) !get memory for PROPP
   mat%propp = 0d0                      !initializes PROPP
   CALL mat_pla(mat%propp,mat%matdef, & !read plastic parameters
                mat%chead,mat%ctail,mat%matdef(12),mat%prope)
   CALL cohill( mat%propp, mat%matdef(3), mat%matdef(8) )  !comput HILL coefficients
 ELSE                            !for elastic material
   mat%matdef(10)  = 1                   !set size to 1 (not used, minimum size)
   ALLOCATE( mat%propp(mat%matdef(10)) )
 END IF

 mat%matdef(11) = 1                      !size of array PROPS (not used, minimum size)
 ALLOCATE( mat%props(mat%matdef(11)) )

 RETURN
 END SUBROUTINE minp_01
