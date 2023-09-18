 SUBROUTINE sc_beam( sec )
 !            --------------------------------
 !            input material properties for
 !            shear deformable 3-D beam (BEAME)
 !            constitutive model

 IMPLICIT NONE
 TYPE (section), POINTER :: sec

 ! local variables
 REAL(kind=8), PARAMETER :: pi = 3.1415926535897932384626433832795
 INTEGER (kind=4) :: matsa(2) = (/ 1,8 /) ! material types allowed for this section
 INTEGER (kind=4) :: matno
 REAL (kind=8) :: h,b,e,g,d,mu,dens
 LOGICAL :: found,plast,isotro,visco

 sec%secty = 42         !3-D beam

 CALL listen('BEAM  ')  !read a line

 matno    = getint('NUMAT ', 1 ,'!Associated Material ..............')
 isotro   = .NOT.exists('NONHOM')
 ! search for associated material
 CALL mat_search(matno,found,sec%mtbas)  !search base material
 sec%mabas = matno                       !material label
 visco = sec%mtbas%mtype == 8
 IF( visco ) mu = sec%mtbas%prope(11)
 IF(.NOT.found) CALL runend ('BEAM  : NOT FOUND ASSOCIATED MATER ')

 IF( .NOT. ANY( sec%mtbas%mtype == matsa) ) &  !check material type
      CALL runend('BEAM  : Invalid material association')

 sec%secdef(1) =   1  !size of integer array
 IF( isotro ) THEN
   sec%secdef(2) =  21  !size of real array
 ELSE
   sec%secdef(2) =  25  !size of real array
   IF( visco ) sec%secdef(2) =  32  !size of real array
 END IF

 ALLOCATE (sec%iprop(sec%secdef(1)), sec%rprop(sec%secdef(2)))
 sec%iprop = 1 !
 IF( isotro ) sec%iprop = 0
 sec%rprop = 0d0 !initializes

 IF( isotro )THEN
   dens = sec%mtbas%prope(5)    !density
   sec%rprop(19) = dens
   e =  sec%mtbas%prope(1)  !young modulus
   IF(exists('SHEARM'))THEN
     g = getrea('SHEARM',0D0,'!Shear Modulus.....................')
   ELSE
     g = sec%mtbas%prope(3) !shear modulus
   END IF
   b = getrea('B     ',0D0,' Width of a rect. cross section ...')
   h = getrea('H     ',0D0,' Height of a rect. cross section ..')
   d = getrea('DIAM  ',0D0,' Diam. of a circu. cross section ..')
   plast = sec%mtbas%matdef(3) > 1
   IF( plast ) THEN
     IF( sec%mtbas%matdef(4) > 2) CALL runend ('Only linear isot. Hard is allowed for beams')
     sec%rprop(21)= sec%mtbas%propp(1)  !K0
     sec%rprop(17)= sec%mtbas%propp(2)  !K'
     sec%rprop(18)=getrea('TOLER ',0D0,' Tolerance for return mapping algor')
     IF(sec%rprop(21) > 0d0) THEN     !plasticity
       sec%rprop(11)=getrea('LAXF  ',0D0,'!Limit axial force ................')
       sec%rprop(12)=getrea('LSHF2 ',0D0,'!Limit shear force in local 2 ax...')
       sec%rprop(13)=getrea('LSHF3 ',0D0,'!Limit shear force in local 3 ax...')
       sec%rprop(14)=getrea('LTM   ',0D0,'!Limit torsional moment ...........')
       sec%rprop(15)=getrea('LBM2  ',0D0,'!Limit bending moment in local 2 ax')
       sec%rprop(16)=getrea('LBM3  ',0D0,'!Limit bending moment in local 3 ax')
     END IF
   END IF
   IF((b <= 0 .OR. h <= 0) .AND. d <=0 )THEN  !User defined values
     sec%rprop(5)=getrea('AREA  ',0D0,'!Axial area .......................')
     sec%rprop(6)=getrea('SHEAR2',0D0,'!Shear area in local 2 direcction .')
     sec%rprop(7)=getrea('SHEAR3',0D0,'!Shear area in local 3 direcction .')
     sec%rprop(8)=getrea('TORINE',0D0,'!Torsional inertia ................')
     sec%rprop(9)=getrea('INERT2',0D0,'!Moment of inertia about local 2 ax')
     sec%rprop(10)=getrea('INERT3',0D0,'!Moment of inertia about local 3 ax')

   ELSE IF( d > 0 ) THEN !standard circular section
     !       elastic constants
     sec%rprop(5) = pi*d**2/4d0
     sec%rprop(6) = sec%rprop(5) !
     sec%rprop(7) = sec%rprop(5)
     sec%rprop(8) = pi*d**4/32d0
     sec%rprop(9) = sec%rprop(8)/2d0
     sec%rprop(10) = sec%rprop(9)
     !       plastic moduli
     IF( plast ) THEN
       sec%rprop(11) = sec%rprop(5)
       sec%rprop(12) = sec%rprop(6)
       sec%rprop(13) = sec%rprop(7)
       sec%rprop(14) = sec%rprop(8)/d
       sec%rprop(15) = sec%rprop(9)/d
       sec%rprop(16) = sec%rprop(10)/d
     END IF
   ELSE  !standard rectangular section
     !       elastic constants
     sec%rprop(5) = b*h
     sec%rprop(6) = b*h*5d0/6d0
     sec%rprop(7) = b*h*5d0/6d0
     sec%rprop(8) = b*h*(b**2+h**2)/12d0
     sec%rprop(9) = b*h**3/12d0
     sec%rprop(10) = h*b**3/12d0
     !       plastic moduli
     IF( plast ) THEN
       sec%rprop(11) = b*h
       sec%rprop(12) = b*h
       sec%rprop(13) = b*h
       IF(h >= b) THEN
         sec%rprop(14) = b**2*(3d0*h-b)/6d0
       ELSE
         sec%rprop(14) = h**2*(3d0*b-h)/6d0
       END IF
       sec%rprop(15) = b*h**2/0.36d1
       sec%rprop(16) = h*b**2/0.36d1
      END IF
   END IF
   sec%rprop(20) = e                                !
   IF(sec%rprop(18) == 0d0) sec%rprop(18) = 0.1d-5  !tolerance for plasticity

   IF (visco )THEN
     sec%rprop(11) = sec%rprop( 5)*mu
     sec%rprop(12) = sec%rprop( 6)*mu
     sec%rprop(13) = sec%rprop( 7)*mu
     sec%rprop(14) = sec%rprop( 8)*mu
     sec%rprop(15) = sec%rprop( 9)*mu
     sec%rprop(16) = sec%rprop(10)*mu
   END IF
   ! compute inertia    matrix
    sec%rprop( 1) = sec%rprop( 5)*dens
    sec%rprop( 2) = sec%rprop( 8)*dens
    sec%rprop( 3) = sec%rprop( 9)*dens
    sec%rprop( 4) = sec%rprop(10)*dens
   ! compute elasticity matrix
   sec%rprop( 5) = sec%rprop( 5)*e
   sec%rprop( 6) = sec%rprop( 6)*g
   sec%rprop( 7) = sec%rprop( 7)*g
   sec%rprop( 8) = sec%rprop( 8)*g
   sec%rprop( 9) = sec%rprop( 9)*e
   sec%rprop(10) = sec%rprop(10)*e

   ! compute yield values
   IF( plast ) THEN
     sec%rprop( 7) = sec%rprop( 7)*sec%rprop(17)
     sec%rprop( 8) = sec%rprop( 8)*sec%rprop(17)/SQRT(3d0)
     sec%rprop( 9) = sec%rprop( 9)*sec%rprop(17)/SQRT(3d0)
     sec%rprop(10) = sec%rprop(10)*sec%rprop(17)/SQRT(3d0)
     sec%rprop(11) = sec%rprop(11)*sec%rprop(17)
     sec%rprop(12) = sec%rprop(12)*sec%rprop(17)
   END IF

 ELSE  !non-homogeneous material
   ! for mass matrix only
   sec%rprop( 1)=getrea('AXIALD',0D0,'!Axial density ....................')
   sec%rprop( 2)=getrea('TORDEN',0D0,'!Torsional density ................')
   sec%rprop( 3)=getrea('INERT2',0D0,'!Moment of inertia about local 2 ax')
   sec%rprop( 4)=getrea('INERT3',0D0,'!Moment of inertia about local 3 ax')
   IF(sec%rprop(2) == 0d0)sec%rprop(2) =sec%rprop(3)+sec%rprop(4) 
   !
   sec%rprop( 5)=getrea('EA    ',0D0,'!Axial stiffness ..................')
   sec%rprop( 6)=getrea('GAY   ',0D0,'!Shear Y Stiffness ................')
   sec%rprop( 7)=getrea('GAZ   ',0D0,'!Shear Z Stiffness ................')
   sec%rprop( 8)=getrea('GJ    ',0D0,'!Torsional stiffness ..............')
   sec%rprop( 9)=getrea('EIY   ',0D0,'!Bending Y Stiffness ..............')
   sec%rprop(10)=getrea('EIZ   ',0D0,'!Bending Z Stiffness...............')
   sec%rprop(11)=getrea('CNTY  ',0D0,' Coupling N-TY.....................')
   sec%rprop(12)=getrea('CNTZ  ',0D0,' Coupling N-TZ.....................')
   sec%rprop(13)=getrea('CNMX  ',0D0,' Coupling N-MX.....................')
   sec%rprop(14)=getrea('CNMY  ',0D0,' Coupling N-MY.....................')
   sec%rprop(15)=getrea('CNMZ  ',0D0,' Coupling N-MZ.....................')
   sec%rprop(16)=getrea('CTYTZ ',0D0,' Coupling TY-TZ....................')
   sec%rprop(17)=getrea('CTYTZ ',0D0,' Coupling TY-TZ....................')
   sec%rprop(18)=getrea('CTYMX ',0D0,' Coupling TY-MX....................')
   sec%rprop(19)=getrea('CTYMY ',0D0,' Coupling TY-MY....................')
   sec%rprop(20)=getrea('CTZMX ',0D0,' Coupling TZ-MX....................')
   sec%rprop(21)=getrea('CTZMY ',0D0,' Coupling TZ-MY....................')
   sec%rprop(22)=getrea('CTZMZ ',0D0,' Coupling TZ-MZ....................')
   sec%rprop(23)=getrea('CMXMY ',0D0,' Coupling MX-MY....................')
   sec%rprop(24)=getrea('CMXMZ ',0D0,' Coupling MX-MZ....................')
   sec%rprop(25)=getrea('CMYMZ ',0D0,' Coupling MY-MZ....................')
   IF (visco )THEN
     e = sec%mtbas%prope(1)  !young modulus
     g = sec%mtbas%prope(3)  !shear modulus
     sec%rprop(26) = sec%rprop( 5)*mu/e
     sec%rprop(27) = sec%rprop( 6)*mu/g
     sec%rprop(28) = sec%rprop( 7)*mu/g
     sec%rprop(29) = sec%rprop( 8)*mu/g
     sec%rprop(30) = sec%rprop( 9)*mu/e
     sec%rprop(31) = sec%rprop(10)*mu/e
     sec%rprop(32) = mu/e
   END IF
 END IF

  RETURN

 END SUBROUTINE sc_beam
