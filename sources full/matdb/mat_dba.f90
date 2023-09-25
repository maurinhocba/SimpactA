 MODULE mat_dba

   ! material data base
   USE param_db,ONLY: mvar, mnam, midn
   USE c_input
   USE flc_db
   IMPLICIT NONE

   INTEGER (kind=4) :: numats = 0! number of materials
   INTEGER (kind=4) :: nusect = 0 ! number of section properties

   ! Derived type for a curve
   TYPE curve
     INTEGER (kind=4) :: np            !number of points
     CHARACTER (len=mnam) :: name        !Associated label
     REAL(kind=8), POINTER :: val(:,:) !x-y values + derivatives
     TYPE(curve), POINTER :: next      !pointer for next curve
   END TYPE curve

   ! Derived type for post-process variables
   TYPE postv
     INTEGER (kind=4) :: type          !type 0:scalar 1:vector 2:tensor
     INTEGER (kind=4) :: dim           !space dimension 1 2 3
     CHARACTER (len=mvar) :: name(7)     !variable name and components
     TYPE(postv), POINTER :: next      !pointer for next variable
   END TYPE postv

   ! Derived type for material
   TYPE mater
     INTEGER (kind=4) :: matno     ! material number (label)
     INTEGER (kind=4) :: mtype     ! material type
                                   ! Presently defined
                                   ! (1) large strain elastic-plastic for metals
                                   !     plastic anisotropy, non-linear isotropic hardening
                                   ! (2)
                                   ! (3)
                                   ! (4) includes thermal properties (element 24 and 19)
                                   ! (5) orthotropic elastic-plastic (space mapping, Oller)
                                   !     under development
                                   ! (6) Hyperelastic & Hyperfoam material
                                   ! (8) Simple visco-elastic material
                                   !(20) Variable material properties, defined by curves
                                   !     under development

     INTEGER (kind=4) :: matdef(12)! (1)Material dimension (1 2 or 3)
                                   ! (2)Elasticity model
                                   !     1: RIGID  Rigid material
                                   !     2: ISOTRO Linear elastic isotropic material
                                   !     3: ORTHOT Linear elastic Orthotropic materials
                                   !     4: OGDEN  Ogden model for rubbers
                                   !     5: POINTS Elastic curve by points
                                   !     6: NONCOM non-compresion
                                   ! (3)Yield criteria
                                   !     1: ELASTI elastic
                                   !     2: MISES  von mises
                                   !     3: HILL48 hill-48
                                   !     4: HILL79 hill-79
                                   !     5: HILL90 hill-90
                                   !     6: BBC03  Banabic 2003 !FR
                                   ! (4)Isotropic hardenig
                                   !     1: IHNONE no Hardening
                                   !     2: IHLINE Linear Hardening
                                   !     3: IHLUDW Ludwik Nadai
                                   !     4: IHSATU linear + saturation law
                                   !     5: IHPOIN defined by points
                                   !     6: IHHOLO Holomon
                                   !     7: IHVOCE Voce (saturation law)
                                   ! (5)Kinematic hardening
                                   !     1: KHNONE no Hardening
                                   !     2: KHLINE linear Hardening
                                   !     3: KHSATU linear + saturation law
                                   !     4: KHPOIN defined by points
                                   !     5: KINDAM
                                   ! (6)Viscous model
                                   !     1: ELASTI elastic
                                   !     2: ELAPLA Elastic-plastic
                                   !     3: ELVIPL Elastic-viscoplastic
                                   ! (7)Thermal model
                                   !     1: ?
                                   ! (8)Sub-model or any other auxiliar value:
                                   !     for MTYPE = 1  (metal)
                                   !       1: ASSOCI Associative plasticity
                                   !       2: NONASO Non-associative plasticity
                                   !     for MTYPE = 2  : submodel
                                   !     for MTYPE = 5 (orthotropic)
                                   !       1234: flag for general orthotropic
                                   !     for MTYPE = 10 (rubber)
                                   !       1-3: number of terms in Ogden model
                                   ! (9)size of array PROPE
                                   !(10)size of array PROPP
                                   !(11)size of array PROPS
                                   !(12)number of curves

     ! values stored in the next arrays depends on the material type (see each material definition)
     REAL (kind=8), POINTER :: prope(:), & !Elasticity properties
                               propp(:), & !Plasticity properties
                               props(:)    !Other properties

     TYPE(curve), POINTER :: chead,ctail   !pointers (first and last) to associated curves

     TYPE (mater), POINTER :: next !pointer to next material in list
   END TYPE mater

   !Derived type for an array of pointers to materials
   TYPE pmater
      TYPE (mater), POINTER :: p
   END TYPE pmater

   TYPE section
     INTEGER (kind=4) :: secno     ! section number (label)
     INTEGER (kind=4) :: secty     ! section type
                                   ! (0) : SOLID (2-D and 3-D solids, TR2D, QUADL, SOLAG, RIGID)
                                   !                                  QUAD4, SOLID
                                   ! (2) : SPOT (point stiffness and damper, SPOT)
                                   ! (5) : SOLSH Layered solid
                                   !(12) : SHELL2 (3-D classic shell,  LBST, CBST)
                                   !(13) : SHELL3 (3-D laminated shell, LBST, CBST)
                                   !(14) : SHELL4 (3-D shear deformable shell, SHELQ, SHELT)
                                   !(41) : TRUSS  (2/3-D truss, TRUSS)
                                   !(42) : BEAM   (3-D shear deformable beam , BEAME)
                                   !(43) : SHREV1 (2-D shear deformable shell/beam , SHREV)
     INTEGER (kind=4) :: mabas     ! material base (label)
     INTEGER (kind=4) :: secdef(5) ! (1) size of array of integer
                                   ! (2) size of array of reals
                                   ! (3) any auxiliar value
                                   ! (4) Number of variables to post-process
                                   ! (5) Number of scalars to post-process
     INTEGER (kind=4),POINTER :: iprop(:) !array of integer
     REAL (kind=8),   POINTER :: rprop(:) !array of real properties
     TYPE (mater), POINTER :: mtbas       !pointer to associated material
     TYPE (postv), POINTER :: postp       !pointer to first postprocess variable
     TYPE (section), POINTER :: next

   END TYPE section

   !Derived type for an array of pointers to sections
   TYPE psection
      TYPE (section), POINTER :: p
   END TYPE psection


   TYPE (mater), POINTER :: mhead, mtail
   TYPE (pmater), POINTER :: pmats(:)    !array of pointers to materials

   TYPE (section), POINTER :: shead, stail
   TYPE (psection), POINTER :: psecs(:)  !array of pointers to sections

   ! scratch for IMPORT facility
   INTEGER(kind=4) :: snn(2,20)


 CONTAINS

   SUBROUTINE ini_mat (mat)

     ! Initializes a material

     TYPE (mater), POINTER :: mat

     ALLOCATE (mat)       !get memory
     NULLIFY (mat%prope, mat%propp, mat%props ) !nullify all real pointers
     mat%matno = 0            !give null or default values to all variables
     mat%mtype = 0
     mat%matdef = (/ 3, 0, 0, 0, 0, 0, 0, 0,&
                     1, & !size of elasticity properties array PROPE
                     1, & !size of plasticity properties array PROPP
                     1, & !size of other properties array PROPS
                     0 /) !number of associated curves
     NULLIFY (mat%chead, mat%ctail, mat%next)
     RETURN
   END SUBROUTINE ini_mat

   SUBROUTINE ini_sect (sect)

     ! Initializes a material

     TYPE (section), POINTER :: sect

     ALLOCATE (sect)       !get memory
     NULLIFY (sect%iprop, sect%rprop ) !nullify all pointers
     sect%secno = 0            !give null or default values to all variables
     sect%secty = 0
     sect%mabas = 0
     sect%secdef = (/ 1, & !size of array of integer
                      1, & !size of array of reals
                      0,0,0 /) !anything, no variables to post-process
     NULLIFY (sect%mtbas, sect%next) !nullify pointer to associated material and next section
     RETURN
   END SUBROUTINE ini_sect

   SUBROUTINE mat_search (matno,found,mat,mid)

     ! searches for material index MATNO (number) in MAT_DB
     !

     LOGICAL, INTENT(OUT) :: found
     INTEGER (kind=4), INTENT(IN) :: matno
     TYPE(mater), POINTER :: mat  !point to material
     INTEGER (kind=4), OPTIONAL :: mid

     INTEGER (kind=4) :: i

     found = .TRUE.          !initializes

     mat => mhead
     i = 0
     DO
       IF (.NOT. ASSOCIATED (mat)) EXIT
       i = i+1
       IF (mat%matno == matno) THEN
         IF( PRESENT(mid) )mid = i
         RETURN
       END IF
       mat => mat%next
     END DO
     found = .FALSE.         !all material  tested and NOT found
     RETURN
   END SUBROUTINE mat_search

   SUBROUTINE sect_search (secno,found,sect,sid)

     ! searches for section index SECNO (number) in SECT_DB
     !

     LOGICAL, INTENT(OUT) :: found
     INTEGER (kind=4), INTENT(IN) :: secno
     TYPE(section), POINTER :: sect
     INTEGER (kind=4), OPTIONAL :: sid

     INTEGER (kind=4) :: i

     found = .TRUE.          !initializes

     sect => shead
     i = 0
     DO
       i = i+1
       IF (.NOT. ASSOCIATED (sect)) EXIT
       IF (sect%secno == secno) THEN
         IF( PRESENT(sid) )sid = i
         RETURN
       END IF
       sect => sect%next
     END DO
     found = .FALSE.         !all sections tested and NOT found

     RETURN
   END SUBROUTINE sect_search

   SUBROUTINE add_cur (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (curve), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add a curve to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     END IF
   END SUBROUTINE add_cur


   SUBROUTINE add_mat (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (mater), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add a material to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     END IF
   END SUBROUTINE add_mat


   SUBROUTINE add_sect (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (section), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add a section to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     END IF
   END SUBROUTINE add_sect

!**************************************************************
   FUNCTION get_nvar (mat,ntype)
     ! get number of internal variables

     IMPLICIT NONE
     INTEGER (kind=4) :: get_nvar
     TYPE(mater), POINTER :: mat
     INTEGER :: ntype

     get_nvar = 0 !elastic material
     IF( mat%matdef(3) <= 1 ) RETURN  !elastic

     SELECT CASE (ntype)
!------------------------------------------------
     CASE (0)        !3D solid
!------------------------------------------------
       SELECT CASE (mat%mtype)
       CASE DEFAULT
         get_nvar = 7
         IF( mat%matdef(5) > 1 ) get_nvar = 13  !kinematic hardening
       END SELECT
!------------------------------------------------
     CASE (1)        !plane stress
!------------------------------------------------
       SELECT CASE (mat%mtype)
       CASE (1,4,5,6)
         get_nvar = 4
       CASE (2,3)
         get_nvar = 4
         IF( mat%matdef(5) > 1 ) get_nvar = 7  !kinematic hardening
       END SELECT
!------------------------------------------------
    !CASE (2:3)      !plane strain & axilsymmetric
!------------------------------------------------
     !  get_nvar = 5
     !CASE (4)
     !  get_nvar = 6
!------------------------------------------------
     END SELECT

   END FUNCTION get_nvar
!**************************************************************


!**************************************************************
 SUBROUTINE elastiff(c,mat,ntype)
 !provides elasticity matrix (vector form)

 REAL (kind=8) :: c(:)
 TYPE (mater), POINTER :: mat
 INTEGER :: ntype
 REAL(kind=8) :: c1,c2,nu,g,ep(9)
 CHARACTER(75) :: warntext
 warntext='WARNING for developer: please check subroutine ELASTIFF!'

 SELECT CASE (ntype)
!------------------------------------------------
 CASE (0)        !3D solid
!------------------------------------------------
   SELECT CASE (mat%mtype)
   CASE (1)
     ! In case of elastic material, C36 is not stored
     ! construct orth properties vector:
     ep(1:3) = mat%prope(1)
     ep(4:6) = mat%prope(2)
     ep(7:9) = mat%prope(3)
     ! Construct stiffness matrix
     CALL omat3d(c,ep,.TRUE.,.FALSE.)
   CASE (5)
     c(1:36) = mat%prope(45:80)
   CASE DEFAULT
     c(1:36) = mat%prope(45:80)
     WRITE(*) warntext
     WRITE(lures,ERR=9999) warntext
   END SELECT
!------------------------------------------------
 CASE (1)        !plane stress  (C11,C12,C22,C33)
   SELECT CASE (mat%mtype)
   CASE (1)
     c(1:4) = mat%prope(7:10)
   CASE (4)
     c(1:4) = mat%prope(10:13)
   CASE (5)
     c(1:4) = mat%prope(16:19)
   CASE (6)
     c(1:4) = mat%prope(19:22)
   CASE DEFAULT
     c(1:4) = mat%prope(16:19)
     WRITE(*) warntext
     WRITE(lures,ERR=9999) warntext
   END SELECT
!------------------------------------------------
 !CASE (2:3)      !plane strain & axilsymmetric
 !  c(1:21) =
!------------------------------------------------
 CASE (4)   !plane stress + shear values
   SELECT CASE (mat%mtype)
   CASE (1)
     c(1:4) = mat%prope(7:10)
     c(5:6) = mat%prope(3)  !*(/ 1d0, 1d0 /)
   !CASE (4)
   !  c(1:4) = mat%prope(10:13)
   CASE (5)
     c(1:4) = mat%prope(16:19)
     c(5:6) =  (/ mat%prope(14),mat%prope(15)  /)
   !CASE (6)
   !  c(1:4) = mat%prope(19:22)
   CASE DEFAULT
     c(1:4) = mat%prope(16:19)
     WRITE(*) warntext
     WRITE(lures,ERR=9999) warntext
   END SELECT
 CASE (5)        !3D laminated solid-shell
!------------------------------------------------
   SELECT CASE (mat%mtype)
   CASE (1)
     ! In case of elastic material, C36 is not stored
     ! construct properties vector:
     nu = mat%prope(2)                        !Poisson ratio
     g  = mat%prope(3)                        !Shear Modulus
     c1 = 2d0*g*(1d0-nu)/(1d0-2d0*nu)        !Diagonal coefficient
     c2 = 2d0*g*nu/(1d0-2d0*nu)              !Off-diagonal Coefficient
     c(1:9) = (/c1,c2,c2,c1,c2,c1,g,g,g/)
   CASE (5)
     c(1:9) = mat%prope((/45,46,47,52,53,59,66,73,80/))
   !CASE DEFAULT
   !  c(1:36) = mat%prope(45:80)
   END SELECT
!------------------------------------------------
 CASE DEFAULT
     WRITE(*) warntext
     WRITE(lures,ERR=9999) warntext
!------------------------------------------------
 END SELECT
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elastiff
!**************************************************************


   SUBROUTINE dump_mate
   !   dumps material database and section database
   IMPLICIT NONE

   INTEGER (kind=4) :: i,j
   TYPE(curve), POINTER :: cur
   TYPE(mater), POINTER :: mat
   TYPE(section), POINTER :: sect
   TYPE(postv), POINTER :: postp
   !------------------------------------------------------------

   WRITE (50,ERR=9999) numats,nusect     !number of materials and sections

   mat => mhead          !point to first material
   DO i = 1, numats      ! for each material
     WRITE (50,ERR=9999) mat%matno, mat%mtype, mat%matdef !write all the control variables
     WRITE (50,ERR=9999) mat%prope, mat%propp, mat%props  !write all the properties
     ! write associated curves if exists
     IF( mat%matdef(12) > 0 )THEN
       cur => mat%chead         !point to first curve
       DO j=1,mat%matdef(12)   !for each curve
         WRITE(50,ERR=9999) cur%np       !number of points in the curve
         WRITE(50,ERR=9999) cur%val      !curve values
         cur => cur%next       !point to next curve
       END DO
     END IF
     mat => mat%next           !point to next material
   END DO

   sect => shead         ! point to first section
   DO i = 1, nusect      ! for each section
     WRITE (50,ERR=9999) sect%secno, sect%secty, sect%mabas, sect%secdef !write all the control variables
     WRITE (50,ERR=9999) sect%iprop, sect%rprop !write all the properties
     IF( sect%secdef(4) > 0)THEN
       postp => sect%postp
       DO j=1,sect%secdef(4)
         WRITE (50,ERR=9999) postp%type,postp%dim,postp%name
         postp => postp%next
       END DO
     END IF
     sect => sect%next
   END DO

   RETURN
   9999 CALL runen2('')
   END SUBROUTINE dump_mate

   SUBROUTINE rest_mate

   !   restores material database

   IMPLICIT NONE

   INTEGER (kind=4) :: i,j
   TYPE(curve), POINTER :: cur
   TYPE(mater), POINTER :: mat
   TYPE(section), POINTER :: sect
   TYPE(postv), POINTER :: postp,post1
   LOGICAL :: found
   !------------------------------------------------------------

   READ (51) numats,nusect     !read number of materials and sections

   NULLIFY (mhead, mtail)

   DO i = 1, numats            !for each material
     ALLOCATE (mat)
     READ (51) mat%matno, mat%mtype, mat%matdef !read all the control variables
     ALLOCATE (mat%prope(mat%matdef( 9)), &
               mat%propp(mat%matdef(10)), &
               mat%props(mat%matdef(11)))
     READ (51) mat%prope, mat%propp, mat%props  !read all the properties
     ! read associated curves if exists
     IF( mat%matdef(12) > 0 )THEN
       NULLIFY (mat%chead, mat%ctail)
       DO j=1,mat%matdef(12)   !for each curve
         ALLOCATE(cur)
         READ(51)cur%np            !number of points in the curve
         ALLOCATE( cur%val(3,cur%np) ) !get memory
         READ(51)cur%val           !curve values
         CALL add_cur(cur,mat%chead,mat%ctail) !add curve to the list
       END DO
     END IF
     CALL add_mat( mat, mhead, mtail )  !add material to the list

   END DO

   ALLOCATE (pmats(numats))
   mat => mhead
   DO i=1,numats
     pmats(i)%p => mat
     mat => mat%next
   END DO

   NULLIFY (shead, stail)

   DO i = 1, nusect            !for each section
     ALLOCATE (sect)
     READ (51) sect%secno, sect%secty, sect%mabas, sect%secdef !read all the control variables
     ALLOCATE (sect%iprop(sect%secdef(1)), &
               sect%rprop(sect%secdef(2)))
     READ (51) sect%iprop, sect%rprop      !read all the properties
     IF( sect%secdef(4) > 0)THEN           !check if list of variables is present
       ALLOCATE (sect%postp)               !get memory for first variable
       postp => sect%postp                 !point auxiliar
       j = 0                               !initializes counter
       DO                                  !loop
         j = j+1                           !increase loop counter
         READ (51) postp%type,postp%dim,postp%name   !read data
         IF( j == sect%secdef(4) )EXIT     !if all variables read, exit
         ALLOCATE (post1)                  !get memory for next variable
         postp%next => post1               !put in the list
         postp => postp%next               !go to next
       END DO
     END IF
     IF( sect%mabas >= 0 )CALL mat_search(sect%mabas,found,sect%mtbas)  !search base material
     CALL add_sect( sect, shead, stail )   !add section to the list

   END DO

   ALLOCATE (psecs(nusect))
   sect => shead
   DO i=1,nusect
     psecs(i)%p => sect
     sect => sect%next
   END DO
   RETURN

   END SUBROUTINE rest_mate

   INCLUDE 'cohill.fi'
   INCLUDE 'conv_smeas.fi'
   INCLUDE 'curaux.fi'
   INCLUDE 'mat_ela.fi'
   INCLUDE 'mat_inp.fi'
   INCLUDE 'mat_pla.fi'
   INCLUDE 'minp_01.fi'
   INCLUDE 'minp_04.fi'
   INCLUDE 'minp_05.fi'
   INCLUDE 'minp_06.fi'
   INCLUDE 'minp_08.fi'
   INCLUDE 'minp_20.fi'
   INCLUDE 'minp_30.fi'
   INCLUDE 'mtype_c.fi'
   INCLUDE 'rdmpcr.fi'
   INCLUDE 'sc_beam.fi'
   INCLUDE 'sc_shel2.fi'
   INCLUDE 'sc_shel3.fi'
   INCLUDE 'sc_shel4.fi'
   INCLUDE 'sc_shel5.fi'
   INCLUDE 'sc_shre1.fi'
   INCLUDE 'sc_solid.fi'
   INCLUDE 'sc_solsh.fi'
   INCLUDE 'sc_spot.fi'
   INCLUDE 'sc_truss.fi'
   INCLUDE 'wrt_mat_db.fi'

 END MODULE mat_dba
