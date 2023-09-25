 SUBROUTINE sc_shel4( sec )
 !            --------------------------------
 !            input material properties for
 !            shear deformable shell (SHELQ, SHELT)
 !            constitutive model

 IMPLICIT NONE
 TYPE (section), POINTER :: sec

 !SECDEF(3) = 0   layered isotropic
 !          = 1   stress resultant isotropic with large strains
 !          = 2   stress resultant isotropic with small strains
 !          = 3   SPECIAL orthotropic
 !          = 4   GENERAL orthotropic
 !          = 5   GENERAL orthotropic with COUPLING

 !IPROP (1) = nlayr   (2) = nothing yet
 !Size of RPROP   12 for elastic isotropic material
 !                37 for plastic isotropic material
 !                15 for elastic special orthotropic material
 !                20 for elastic symmetric orthotropic material
 !                41 for elastic general orthotropic material

 !RPROP(1:5) = thick,min_thick,max_thick, ?,s_dens

 !local variables
 INTEGER (kind=4), PARAMETER :: matsa(2) = (/ 1,5 /) ! material types allowed for this section
 INTEGER (kind=4) :: matno,nlayr,rtype,i
 REAL (kind=8) :: thick,auxil,young,poiss,f,g,h,dens
 LOGICAL :: found,readm,plast
 TYPE (mater), POINTER :: mat

 sec%secty = 14         !shell type 4

 CALL listen('SHELL4')  !read a line

 ! associate a material
 matno    = getint('NUMAT ', 1 ,'!Associated Material ..............')
 CALL mat_search(matno,found,sec%mtbas)  !search base material
 sec%mabas = matno                       !material label
 IF(.NOT.found) CALL runend ('SHELL4: NOT FOUND ASSOCIATED MATER ')
 mat => sec%mtbas

 IF( .NOT. ANY( mat%mtype == matsa) ) &  !check if associated material is valid
      CALL runend('SHELL4: Invalid material association')

 readm = exists('READMA')  !read elastic matrix ?
 plast = .FALSE.           !initializes to Elastic

 nlayr = getint('NLAYR ',   1,' NUMBER OF LAYERS .................')
 IF ( nlayr  <= 0) CALL runend('SHELL4: VERIFY NUMBER OF LAYERS ')
 thick = getrea('THICK ',1.D0,'!Thickness ........................')
 IF ( thick <= 0.D0) CALL runend('SHELL4: VERIFY THICKNESS ')
 dens = mat%prope(5)

 sec%secdef(1) =  2  !size of integer array
 SELECT CASE ( mat%mtype )  !according to material base type
 CASE (1)                   !metals
   sec%secdef(2) =  12            !size of real array
   IF( nlayr == 1 )THEN           !stress resultant model
     IF( exists('LARGE ') )THEN
       sec%secdef(3) =  1     !stress resultant isotropic with large strains
     ELSE
       sec%secdef(3) =  2     !stress resultant isotropic with small strains
     END IF
   ELSE
     sec%secdef(3) =  0       !layered isotropic elastic, orthotropic plastic
   END IF
   IF( mat%matdef(3) > 1 ) THEN  !base material is plastic
     plast = .TRUE.                    !modify flag
     sec%secdef(2) =  37               !size of real array including plastic parameters
   END IF
   IF( readm ) THEN
     rtype = 1   !isotropic with user defined matrix
   ELSE
     rtype = 0   !isotropic with standard elastic matrix
   END IF

 CASE (5)                   !Orthotropic material
   ! default values
   rtype = 2                     !SPECIAL orthotropic standard
   sec%secdef(2) =  15           !size of real array
   sec%secdef(3) =  3            !SPECIAL orthotropic
   IF( readm )THEN
     IF( exists('GENERA') )THEN  !general orthotropic matrix
       sec%secdef(2) =  20       !size of real array
       sec%secdef(3) =  4        !GENERAL orthotropic
       rtype = 4              !read general matrix
     ELSE IF( exists('COUPLI') )THEN      !coupling included
       sec%secdef(2) =  41       !size of real array
       sec%secdef(3) =  5        !GENERAL orthotropic with COUPLING
       rtype = 5                 !read orthotropic with coupling
     ELSE
       rtype = 3              !read SPECIAL orthotropic
     END IF
   END IF
 CASE (9)                     !not yet
   sec%secdef(2) =  1
 END SELECT

 ALLOCATE (sec%iprop(sec%secdef(1)), sec%rprop(sec%secdef(2)))
 !  read parameters or compute elastic matrix
 IF( nlayr > 1 .AND. rtype > 0 )CALL runend &    !check consistency
    ('SHEL4: nlayr > 1 only valid for isotropic materials')
 IF( readm )CALL listen('SHEL4 ')
 SELECT CASE (rtype)
 CASE (0)  !does not read a matrix
   IF( nlayr == 1 )THEN     ! integrated constitutive matrix
     sec%rprop( 6) = sec%mtbas%prope( 7)*thick
     sec%rprop( 7) = sec%mtbas%prope( 8)*thick
     sec%rprop( 8) = sec%mtbas%prope(10)*thick
     auxil = thick*thick/12d0
     sec%rprop( 9) = sec%rprop( 6)*auxil
     sec%rprop(10) = sec%rprop( 7)*auxil
     sec%rprop(11) = sec%rprop( 8)*auxil
     sec%rprop(12) = 5d0*sec%rprop( 8)/6d0
   ELSE                     !layered
     sec%rprop(6:9) = sec%mtbas%prope(7:10)   !plane stress elasticity matrix
     sec%rprop(10 ) = sec%rprop(9)*5d0/6d0    !transverse shear modulus
     sec%rprop(11 ) = thick                   !thickness
     sec%rprop(12 ) = thick*thick/12d0        !factor  !OVERLAYED
   END IF
 CASE (1)   !isotropic laminate
 !             6  7
 !                6
 !                   8
 !                      9 10
 !                         9
 !                           11
 !                              12
 !                                 12
   sec%rprop( 6)=getrea('C1111 ',0d0,'!Membrane Coeff dir 1 .............')
   sec%rprop( 7)=getrea('C1122 ',0d0,'!Membrane Coeff dir 1-2............')
   sec%rprop( 8)=getrea('C1212 ',0d0,'!Membrane Shear Coeff 1-2..........')
   sec%rprop( 9)=getrea('D1111 ',0d0,'!Bending  Coeff dir 1 .............')
   sec%rprop(10)=getrea('D1122 ',0d0,'!Bending  Coeff dir 1-2............')
   sec%rprop(11)=getrea('D1212 ',0d0,'!Bending  Shear Coeff 1-2..........')
   sec%rprop(12)=getrea('G1313 ',0d0,'!Transver Shear Coeff 1-3..........')
 CASE (2)   !standard special orthotropic laminate (does not read a matrix)
 !             6  7
 !                8
 !                   9
 !                     10 11
 !                        12
 !                           13
 !                              14
 !                                 15
   sec%rprop( 6)=mat%prope(20)*thick           !Membrane Coeff dir 1 .....
   sec%rprop( 7)=mat%prope(21)*thick           !Membrane Coeff dir 1-2....
   sec%rprop( 8)=mat%prope(24)*thick           !Membrane Coeff dir 2 .....
   sec%rprop( 9)=mat%prope(28)*thick           !Membrane Shear Coeff 1-2..
   auxil = thick*thick/12d0
   sec%rprop(10)=sec%rprop( 6)*auxil           !Bending  Coeff dir 1 .....
   sec%rprop(11)=sec%rprop( 7)*auxil           !Bending  Coeff dir 1-2....
   sec%rprop(12)=sec%rprop( 8)*auxil           !Bending  Coeff dir 2 .....
   sec%rprop(13)=sec%rprop( 9)*auxil           !Bending  Shear Coeff 1-2..
   sec%rprop(14)=mat%prope(73)*thick*5d0/6d0   !Transver Shear Coeff 1-3..
   sec%rprop(15)=mat%prope(80)*thick*5d0/6d0   !Transver Shear Coeff 2-3..
 CASE (3)   !Special orthotropic laminate
   sec%rprop( 6)=getrea('C1111 ',0d0,'!Membrane Coeff dir 1 .............')
   sec%rprop( 7)=getrea('C1122 ',0d0,'!Membrane Coeff dir 1-2............')
   sec%rprop( 8)=getrea('C2222 ',0d0,'!Membrane Coeff dir 2 .............')
   sec%rprop( 9)=getrea('C1212 ',0d0,'!Membrane Shear Coeff 1-2..........')
   sec%rprop(10)=getrea('D1111 ',0d0,'!Bending  Coeff dir 1 .............')
   sec%rprop(11)=getrea('D1122 ',0d0,'!Bending  Coeff dir 1-2............')
   sec%rprop(12)=getrea('D2222 ',0d0,'!Bending  Coeff dir 2 .............')
   sec%rprop(13)=getrea('D1212 ',0d0,'!Bending  Shear Coeff 1-2..........')
   sec%rprop(14)=getrea('G1313 ',0d0,'!Transver Shear Coeff 1-3..........')
   sec%rprop(15)=getrea('G2323 ',0d0,'!Transver Shear Coeff 2-3..........')
 CASE (4,5)  !General orthotropic laminate
 !             6  7  8
 !                9 10
 !                  11
 !                     12 13 14
 !                        15 16
 !                           17
 !                              18 19
 !                                 20
   sec%rprop( 6)=getrea('C1111 ',0d0,'!Membrane Coeff dir 1 .............')
   sec%rprop( 7)=getrea('C1122 ',0d0,'!Membrane Coeff dir 1-2............')
   sec%rprop( 8)=getrea('C1112 ',0d0,' Membrane Coupling  11-12..........')
   sec%rprop( 9)=getrea('C2222 ',0d0,'!Membrane Coeff dir 2 .............')
   sec%rprop(10)=getrea('C2212 ',0d0,' Membrane Coupling  22-12..........')
   sec%rprop(11)=getrea('C1212 ',0d0,'!Membrane Shear Coeff 1-2..........')
   sec%rprop(12)=getrea('D1111 ',0d0,'!Bending  Coeff dir 1 .............')
   sec%rprop(13)=getrea('D1122 ',0d0,'!Bending  Coeff dir 1-2............')
   sec%rprop(14)=getrea('D1112 ',0d0,' Bending Coupling  11-12 ..........')
   sec%rprop(15)=getrea('D2222 ',0d0,'!Bending  Coeff dir 2 .............')
   sec%rprop(16)=getrea('D2212 ',0d0,' Bending Coupling  22-12 ..........')
   sec%rprop(17)=getrea('D1212 ',0d0,'!Bending  Shear Coeff 1-2..........')
   sec%rprop(18)=getrea('G1313 ',0d0,'!Transver Shear Coeff 1-3..........')
   sec%rprop(19)=getrea('G1323 ',0d0,' Transver Shear Coupling 13-23 ....')
   sec%rprop(20)=getrea('G2323 ',0d0,'!Transver Shear Coeff 2-3..........')
   IF( rtype == 5 )THEN
     WRITE (lures,"('  Coupling terms ')",ERR=9999)
     sec%rprop(21)=getrea('B1111 ',0d0,' Coupling term N11-K11.............')
     sec%rprop(22)=getrea('B1122 ',0d0,' Coupling term N11-K22.............')
     sec%rprop(23)=getrea('B1112 ',0d0,' Coupling term N11-K12.............')
     sec%rprop(24)=getrea('B2211 ',0d0,' Coupling term N22-K11.............')
     sec%rprop(25)=getrea('B2222 ',0d0,' Coupling term N22-K22.............')
     sec%rprop(26)=getrea('B2212 ',0d0,' Coupling term N22-K12.............')
     sec%rprop(27)=getrea('B1211 ',0d0,' Coupling term N12-K11.............')
     sec%rprop(28)=getrea('B1222 ',0d0,' Coupling term N12-K22.............')
     sec%rprop(29)=getrea('B1212 ',0d0,' Coupling term N12-K12.............')
     !sec%rprop(30)=getrea('B1113 ',0d0,' Coupling term N11-E13.............')
     !sec%rprop(31)=getrea('B1123 ',0d0,' Coupling term N11-E23.............')
     !sec%rprop(32)=getrea('B2213 ',0d0,' Coupling term N22-E13.............')
     !sec%rprop(33)=getrea('B2223 ',0d0,' Coupling term N22-E23.............')
     !sec%rprop(34)=getrea('B1213 ',0d0,' Coupling term N12-E13.............')
     !sec%rprop(35)=getrea('B1223 ',0d0,' Coupling term N12-E23.............')
     !sec%rprop(36)=getrea('A1113 ',0d0,' Coupling term M11-E13.............')
     !sec%rprop(37)=getrea('A1123 ',0d0,' Coupling term M11-E23.............')
     !sec%rprop(38)=getrea('A2213 ',0d0,' Coupling term M22-E13.............')
     !sec%rprop(39)=getrea('A2223 ',0d0,' Coupling term M22-E23.............')
     !sec%rprop(40)=getrea('A1213 ',0d0,' Coupling term M12-E13.............')
     !sec%rprop(41)=getrea('A1223 ',0d0,' Coupling term M12-E23.............')
   END IF
 END SELECT
 sec%iprop(1) = nlayr      !number of layers
 sec%iprop(2) = 1      !nothing yet

 sec%rprop(1) = thick
 sec%rprop(2) = getrea('MINTHI',0.D0,' Mimimum Thickness Ratio...........')
 sec%rprop(3) = getrea('MAXTHI',1.D3,' Maximum Thickness Ratio...........')
 sec%rprop(4) = 0d0  !nothing
 IF( readm .AND. exists('DENSI',i) )THEN
   sec%rprop(5) = param(i)
 ELSE
   sec%rprop(5) = dens*thick
 END IF
 ! plastic parameters for stress resultant solid isotropic
 ! variables for plasticity
 IF( .NOT. plast ) RETURN
 ! next only for isotropic materials
 IF( mat%matdef(4) > 4 )CALL runend  ('SHELL4: Invalid isotropic hardening')
 IF( nlayr == 1 )THEN                 !analitically Integrated
   ! coeficients for radial RETURN
   young = mat%prope(1)
   poiss = mat%prope(2)
   sec%rprop(21) = mat%propp(1) !conso
   sec%rprop(22) = mat%propp(2) !efref
   sec%rprop(23) = mat%propp(3) !consn
   sec%rprop(24) = mat%propp(4) !refls
   sec%rprop(25) = mat%propp(6) !kinet

   sec%rprop(25) = sec%rprop(24)*thick*2d0/3d0               !khard
   sec%rprop(26) = thick*thick/12d0                          !f1
   sec%rprop(27) = (1d0/thick)**2                            !coefn
   sec%rprop(28) = (4d0/thick/thick)**2                      !coefm
   sec%rprop(29) = 2d0/SQRT(3d0)/thick**3                    !coenm
   sec%rprop(30) = young*thick/(1-poiss)                     !ddmat1
   sec%rprop(31) = young*thick/(1+poiss)                     !ddmat2
   sec%rprop(32) = young*thick/(1+poiss)/2.0d0               !ddmat3
   sec%rprop(33) = 5.0d0/6.0d0*sec%rprop(32)                 !gkh
   sec%rprop(34) = (sec%rprop(30)+sec%rprop(25))/3d0         !pdmat1
   sec%rprop(35) = (sec%rprop(31)+sec%rprop(25))*1d0         !pdmat2
   sec%rprop(36) = (sec%rprop(32)+sec%rprop(25))*2d0         !pdmat3
   sec%rprop(37) = (sec%rprop(33)+sec%rprop(25))*2d0/thick**2  !zeta

 ELSE                                !integrated by layers
   sec%rprop(12) = mat%propp(1) !conso
   sec%rprop(13) = mat%propp(2) !efref
   sec%rprop(14) = mat%propp(3) !consn
   sec%rprop(15) = mat%propp(4) !refls
   sec%rprop(16) = mat%propp(6) !kinet
   !           compute B (flow rull) matrix from Hill coeficients
   f = mat%propp(22)
   g = mat%propp(23)
   h = mat%propp(24)
   sec%rprop(23) = g + h
   sec%rprop(24) =   - h
   sec%rprop(25) = f + h
   sec%rprop(26) = 2d0*mat%propp(25)       !2n
   sec%rprop(27) = 2d0*mat%propp(26)       !2l
   sec%rprop(28) = 2d0*mat%propp(27)       !2m
   !           compute A (yield) matrix from Hill coeficients
   IF( exists('ASSOCIA') ) THEN
     sec%rprop(17:22) = sec%rprop(23:28)
   ELSE
     f = mat%propp(16)
     g = mat%propp(17)
     h = mat%propp(18)
     sec%rprop(17) = g + h
     sec%rprop(18) =   - h
     sec%rprop(19) = f + h
     sec%rprop(20) =  2d0*mat%propp(19)       !2n
     sec%rprop(21) =  2d0*mat%propp(20)       !2l
     sec%rprop(22) =  2d0*mat%propp(21)       !2m
   END IF
 END IF
 sec%rprop(5) = thick*mat%prope(5)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE sc_shel4
