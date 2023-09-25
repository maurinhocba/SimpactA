 SUBROUTINE sc_shre1( sec )
 !            --------------------------------
 !            input material properties for
 !            1 shear deformable 2-D shell/beam (SHREV)
 !            2 classic linear 2-D beam (SHREV NTYPE=4)
 !            3 classic rotation-free 2-D shell/beam (BEAM2)
 !

 IMPLICIT NONE
 TYPE (section), POINTER :: sec

 !SECDEF(3) = 0   layered isotropic

 ! material types allowed for this section
 INTEGER (kind=4) :: matsa(3) = (/ 1,5,6 /)
 INTEGER (kind=4) :: matno,rtype,nlayr
 REAL (kind=8) :: thick,auxil,young,poiss,width,area,inert
 LOGICAL :: found,plast
 TYPE (mater), POINTER :: mat


 sec%secty = 43         !shell of revolution

 CALL listen('SHREV1')  !read a line

 matno    = getint('NUMAT ', 1 ,'!Associated Material ..............')
 ! search for associated material
 CALL mat_search(matno,found,sec%mtbas)  !search base material
 sec%mabas = matno                       !material label
 IF(.NOT.found) CALL runend ('SHREV1: NOT FOUND ASSOCIATED MATER ')

 IF( .NOT. ANY( sec%mtbas%mtype == matsa) ) & !check if the material is correct
      CALL runend('SHREV1: Invalid material association')

 !readm = exists('READMA')  !read elastic matrix ?
 plast = .FALSE.

 nlayr = getint('NLAYR ',   1,' NUMBER OF LAYERS .................')
 IF ( nlayr  <= 0) CALL runend('SHREV1: VERIFY NUMBER OF LAYERS ')
 thick = getrea('THICK ',1.D0,' Thickness ........................')
 IF ( thick == 0.D0) CALL runend('SHREV1: VERIFY THICKNESS ')

 sec%secdef(1) =  2   !size of integer array
 sec%secdef(2) =  35  !size of real array
 SELECT CASE ( sec%mtbas%mtype )
 CASE (1)
   !IF( nlayr == 1 )THEN
     IF( exists('LARGE ') )THEN
       sec%secdef(3) =  1   !stress resultant isotropic with large strains
     ELSE
       sec%secdef(3) =  2   !stress resultant isotropic with small strains
     END IF
   !ELSE
   !  sec%secdef(3) =  0   !layered isotropic
   !END IF
   IF( sec%mtbas%matdef(3) > 1 ) THEN
     plast = .TRUE.
   END IF
   !IF( readm ) THEN
   !  rtype = 1   !isotropic with used defined matrix
   !ELSE
     rtype = 0   !isotropic with standard elastic matrix
   !END IF

 !CASE (5)
 !  ! default values
 !  rtype = 2               !SPECIAL orthotropic standard
 !  sec%secdef(2) =  14     !size of real array
 !  sec%secdef(3) =  3      !SPECIAL orthotropic
 !  IF( readm )THEN
 !    IF( exists('COUPLI') )THEN      !coupling included
 !      sec%secdef(2) =  40 !size of real array
 !      sec%secdef(3) =  4  !GENERAL orthotropic with COUPLING
 !      rtype = 4           !read orthotropic with coupling
 !    ELSE
 !      rtype = 3           !read SPECIAL orthotropic
 !    END IF
 !  END IF
 END SELECT

 ALLOCATE (sec%iprop(sec%secdef(1)), sec%rprop(sec%secdef(2)))

 !  read parameters or compute elastic matrix
 !IF( nlayr > 1 .AND. rtype > 0 )CALL runend &
 !   ('SHEL4: nlayr > 1 only valid for isotropic mats')
 !IF( readm )CALL listen('SHEL4 ')
 SELECT CASE (rtype)
 CASE (0)
   !IF( nlayr == 1 )THEN     !
     young   = sec%mtbas%prope(1)
     poiss   = sec%mtbas%prope(2)

     ! For PLANE STRESS

     !         elastic constants
     width = getrea('WIDTH ',1d0,' Width for 2-D Beam Element........')
     sec%rprop(10) = width
     area  = width*thick
     inert = area*thick**2/12d0
     area = getrea('AREA  ',area,' Section for 2-D Beam Element........')
     inert= getrea('INERT ',inert,' Inertia for 2-D Beam Element........')

     !         elastic matrix coeficients
     auxil = young/(1d0+poiss)/2d0
     sec%rprop(11) = area*young
     sec%rprop(12) = inert*young
     sec%rprop(13) = 5d0/12d0*young*area/(1d0+poiss)

     ! For PLANE STRAIN or AXILSYMMETRIC

     !         elastic matrix coeficients
     auxil = young*thick/(1.0-poiss*poiss)
     sec%rprop(21) = auxil
     sec%rprop(22) = auxil*poiss
     auxil = thick*thick/12d0
     sec%rprop(23) = sec%rprop(21)*auxil
     sec%rprop(24) = sec%rprop(22)*auxil
     sec%rprop(25) = 5d0*thick*young/(1d0+poiss)/12d0
   !ELSE                     !layered
   !  sec%rprop(5:8) = sec%mtbas%prope(7:10)
   !  sec%rprop( 9 ) = sec%rprop(8)*5d0/6d0
   !END IF
 !CASE (1)   !isotropic laminate
 !  sec%rprop( 5)=getrea('C1111 ',0d0,'!Membrane Coeff dir 1 .............')
 !  sec%rprop( 6)=getrea('C1122 ',0d0,'!Membrane Coeff dir 1-2............')
 !  sec%rprop( 7)=getrea('C1212 ',0d0,'!Membrane Shear Coeff 1-2..........')
 !  sec%rprop( 8)=getrea('D1111 ',0d0,'!Bending  Coeff dir 1 .............')
 !  sec%rprop( 9)=getrea('D1122 ',0d0,'!Bending  Coeff dir 1-2............')
 !  sec%rprop(10)=getrea('D1212 ',0d0,'!Bending  Shear Coeff 1-2..........')
 !  sec%rprop(11)=getrea('G1313 ',0d0,'!Transver Shear Coeff 1-3..........')
 !CASE (2)   !standard special orthotropic laminate
 !  sec%rprop( 5)=sec%mtbas%props(1)*thick      !Membrane Coeff dir 1 .....
 !  sec%rprop( 6)=sec%mtbas%props(2)*thick      !Membrane Coeff dir 1-2....
 !  sec%rprop( 7)=sec%mtbas%props(5)*thick      !Membrane Coeff dir 2 .....
 !  sec%rprop( 8)=sec%mtbas%props(9)*thick      !Membrane Shear Coeff 1-2..
 !  auxil = thick*thick/12d0
 !  sec%rprop( 9)=sec%rprop( 5)*auxil           !Bending  Coeff dir 1 .....
 !  sec%rprop(10)=sec%rprop( 6)*auxil           !Bending  Coeff dir 1-2....
 !  sec%rprop(11)=sec%rprop( 7)*auxil           !Bending  Coeff dir 2 .....
 !  sec%rprop(12)=sec%rprop( 8)*auxil           !Bending  Shear Coeff 1-2..
 !  sec%rprop(13)=sec%mtbas%props(54)*thick     !Transver Shear Coeff 1-3..
 !  sec%rprop(14)=sec%mtbas%props(61)*thick     !Transver Shear Coeff 2-3..
 !CASE (3)   !Special orthotropic laminate
 !  sec%rprop( 5)=getrea('C1111 ',0d0,'!Membrane Coeff dir 1 .............')
 !  sec%rprop( 6)=getrea('C1122 ',0d0,'!Membrane Coeff dir 1-2............')
 !  sec%rprop( 7)=getrea('C2222 ',0d0,'!Membrane Coeff dir 2 .............')
 !  sec%rprop( 8)=getrea('C1212 ',0d0,'!Membrane Shear Coeff 1-2..........')
 !  sec%rprop( 9)=getrea('D1111 ',0d0,'!Bending  Coeff dir 1 .............')
 !  sec%rprop(10)=getrea('D1122 ',0d0,'!Bending  Coeff dir 1-2............')
 !  sec%rprop(11)=getrea('D2222 ',0d0,'!Bending  Coeff dir 2 .............')
 !  sec%rprop(12)=getrea('D1212 ',0d0,'!Bending  Shear Coeff 1-2..........')
 !  sec%rprop(13)=getrea('G1313 ',0d0,'!Transver Shear Coeff 1-3..........')
 !  sec%rprop(14)=getrea('G2323 ',0d0,'!Transver Shear Coeff 2-3..........')
 !CASE (4)
 !  sec%rprop( 5)=getrea('C1111 ',0d0,'!Membrane Coeff dir 1 .............')
 !  sec%rprop( 6)=getrea('C1122 ',0d0,'!Membrane Coeff dir 1-2............')
 !  sec%rprop( 7)=getrea('C1112 ',0d0,' Membrane Coupling  11-12..........')
 !  sec%rprop( 8)=getrea('C2222 ',0d0,'!Membrane Coeff dir 2 .............')
 !  sec%rprop( 9)=getrea('C2212 ',0d0,' Membrane Coupling  22-12..........')
 !  sec%rprop(10)=getrea('C1212 ',0d0,'!Membrane Shear Coeff 1-2..........')
 !  sec%rprop(11)=getrea('D1111 ',0d0,'!Bending  Coeff dir 1 .............')
 !  sec%rprop(12)=getrea('D1122 ',0d0,'!Bending  Coeff dir 1-2............')
 !  sec%rprop(13)=getrea('D1112 ',0d0,' Bending Coupling  11-12 ..........')
 !  sec%rprop(14)=getrea('D2222 ',0d0,'!Bending  Coeff dir 2 .............')
 !  sec%rprop(15)=getrea('D2212 ',0d0,' Bending Coupling  22-12 ..........')
 !  sec%rprop(16)=getrea('D1212 ',0d0,'!Bending  Shear Coeff 1-2..........')
 !  sec%rprop(17)=getrea('G1313 ',0d0,'!Transver Shear Coeff 1-3..........')
 !  sec%rprop(18)=getrea('G1323 ',0d0,' Transver Shear Coupling 13-23 ....')
 !  sec%rprop(19)=getrea('G2323 ',0d0,'!Transver Shear Coeff 2-3..........')
 !  WRITE (lures,"('  Coupling terms ')")
 !  sec%rprop(20)=getrea('B1111 ',0d0,' Coupling term N11-K11.............')
 !  sec%rprop(21)=getrea('B1122 ',0d0,' Coupling term N11-K22.............')
 !  sec%rprop(22)=getrea('B1112 ',0d0,' Coupling term N11-K12.............')
 !  sec%rprop(23)=getrea('B2211 ',0d0,' Coupling term N22-K11.............')
 !  sec%rprop(24)=getrea('B2222 ',0d0,' Coupling term N22-K22.............')
 !  sec%rprop(25)=getrea('B2212 ',0d0,' Coupling term N22-K12.............')
 !  sec%rprop(26)=getrea('B1211 ',0d0,' Coupling term N12-K11.............')
 !  sec%rprop(27)=getrea('B1222 ',0d0,' Coupling term N12-K22.............')
 !  sec%rprop(28)=getrea('B1212 ',0d0,' Coupling term N12-K12.............')
 !  sec%rprop(29)=getrea('B1113 ',0d0,' Coupling term N11-E13.............')
 !  sec%rprop(30)=getrea('B1123 ',0d0,' Coupling term N11-E23.............')
 !  sec%rprop(31)=getrea('B2213 ',0d0,' Coupling term N22-E13.............')
 !  sec%rprop(32)=getrea('B2223 ',0d0,' Coupling term N22-E23.............')
 !  sec%rprop(33)=getrea('B1213 ',0d0,' Coupling term N12-E13.............')
 !  sec%rprop(34)=getrea('B1223 ',0d0,' Coupling term N12-E23.............')
 !  sec%rprop(35)=getrea('A1113 ',0d0,' Coupling term M11-E13.............')
 !  sec%rprop(36)=getrea('A1123 ',0d0,' Coupling term M11-E23.............')
 !  sec%rprop(37)=getrea('A2213 ',0d0,' Coupling term M22-E13.............')
 !  sec%rprop(38)=getrea('A2223 ',0d0,' Coupling term M22-E23.............')
 !  sec%rprop(39)=getrea('A1213 ',0d0,' Coupling term M12-E13.............')
 !  sec%rprop(40)=getrea('A1223 ',0d0,' Coupling term M12-E23.............')
 END SELECT
 sec%iprop(1) = nlayr
 sec%iprop(2) = 1

 sec%rprop(1) = thick
 sec%rprop(2) = getrea('MINTHI',0.D0,' Mimimum Thickness Ratio...........')
 sec%rprop(3) = getrea('MAXTHI',1.D3,' Maximum Thickness Ratio...........')
 sec%rprop(4) = area   !area

 ! plastic parameters for stress resultant solid isotropic
 ! variables for plasticity
 IF( .NOT. plast ) RETURN

 mat => sec%mtbas
 !IF (mat%mtype == 1) THEN
   sec%rprop(5:8) = mat%propp(1:4)
   sec%rprop(9) = mat%propp(6)
 !ELSE
 !
 !END IF
 !IF( nlayr == 1 )THEN                     !solid isotropic
   ! coeficients for radial RETURN

 ! For PLANE STRESS

   sec%rprop(14) = thick*thick/12d0
   sec%rprop(15) = 1d0/area**2                                    !coefn
   sec%rprop(16) = 1d0/(area*thick/4d0)**2                        !coefm
   sec%rprop(17) = 2d0/SQRT(3d0)/area**2/thick                    !coenm
   sec%rprop(18) = (young*area+sec%rprop(9))*2d0/3d0              !pdmat1
   sec%rprop(19) = (young*area*5d0/12d0/(1d0+poiss)+sec%rprop(9))*2d0/area**2  !zeta

 ! For PLANE STRAIN or AXILSYMMETRIC

   sec%rprop(26) = thick*thick/12d0                          !f1
   sec%rprop(27) = (1d0/thick)**2                            !coefn
   sec%rprop(28) = (4d0/thick/thick)**2                      !coefm
   sec%rprop(29) = 2d0/SQRT(3d0)/thick**3                    !coenm
   sec%rprop(30) = young*thick/(1-poiss)                     !ddmat1
   sec%rprop(31) = young*thick/(1+poiss)                     !ddmat2
   sec%rprop(32) = 5.0d0/6.0d0*sec%rprop(30)/2d0             !gkh
   sec%rprop(33) = (sec%rprop(30)+sec%rprop(9))/3d0          !pdmat1
   sec%rprop(34) = (sec%rprop(31)+sec%rprop(9))*1d0          !pdmat2
   sec%rprop(35) = (sec%rprop(32)+sec%rprop(9))*2d0/thick**2 !zeta

 !ELSE
 !  sec%rprop(11) = mat%propp(1) !conso
 !  sec%rprop(12) = mat%propp(2) !efref
 !  sec%rprop(13) = mat%propp(3) !consn
 !  sec%rprop(14) = mat%propp(4) !refls
 !  sec%rprop(15) = mat%propp(6) !kinet
 !  !           compute B (flow rull) matrix from Hill coeficients
 !  f = mat%propp(22)
 !  g = mat%propp(23)
 !  h = mat%propp(24)
 !  sec%rprop(22) = g + h
 !  sec%rprop(23) =   - h
 !  sec%rprop(24) = f + h
 !  sec%rprop(25) = 2d0*mat%propp(25)       !2n
 !  sec%rprop(26) = 2d0*mat%propp(26)       !2l
 !  sec%rprop(27) = 2d0*mat%propp(27)       !2m
 !  !           compute A (yield) matrix from Hill coeficients
 !  IF( exists('ASSOCIA') ) THEN
 !    sec%rprop(16:21) = sec%rprop(22:27)
 !  ELSE
 !    f = mat%propp(16)
 !    g = mat%propp(17)
 !    h = mat%propp(18)
 !    sec%rprop(16) = g + h
 !    sec%rprop(17) =   - h
 !    sec%rprop(18) = f + h
 !    sec%rprop(19) =  2d0*mat%propp(19)       !2n
 !    sec%rprop(20) =  2d0*mat%propp(20)       !2l
 !    sec%rprop(21) =  2d0*mat%propp(21)       !2m
 !  END IF
 !END IF

 RETURN
 END SUBROUTINE sc_shre1
