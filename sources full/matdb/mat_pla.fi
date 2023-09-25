 SUBROUTINE mat_pla (ppr, def, chead, ctail, nc, epr)
 !
 !  read standard plasticity constants
 !
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN OUT) :: def(:) ! (3)Yield criteria
                                            ! (4)Isotropic hardenig
                                            ! (5)Kinematic hardening
                                            ! (6)Viscosity
                                            ! (7) ?

 REAL (kind=8), INTENT(OUT) :: ppr(:)   ! properties array
    !(1..4) isotropic hardening parameters
    !(5) death time
    !(6..7) kinematic hardening parameters
    !(8..NN) other properties
    !(16) viscosity (for viscoplastic models)
 TYPE(CURVE), POINTER, OPTIONAL :: chead,ctail
 INTEGER(kind=4), OPTIONAL :: nc
 REAL (kind=8), OPTIONAL :: epr(:)    ! elastic properties array
 ! local variables
 REAL (kind=8) :: uniax,eunif,strem,young !,aux
 LOGICAL :: pars(4),rdcur
 INTEGER (kind=4) :: npar,i
 CHARACTER (len=midn) :: straim, stresm

 IF( def(3) == 1 )RETURN  !elastic

 !IF( def(4) == 5 )RETURN  !IHPOIN, discuss placement of Call rdmpcr

 CALL listen('MATPLA')
 IF( exists('ENDMAT') .OR. exists('SECTIO') .OR. exists('DEFINI'))THEN
   backs = .TRUE.
   def(3) = 1  !elastic
   RETURN
 END IF

 ! first read isotropic hardening law

 rdcur = .FALSE.  !initializes IHPOINT

 IF( def(4) == 0 .AND. def(3) > 1 ) CALL runend('MATPLA: yield data missing ')

 SELECT CASE (def(4))

 CASE (1)  ! IHNONE: No hardening
   ppr(1) = getrea('UNIAX ',0d0,' Yield Stress .....................')
 CASE (2)  ! IHLINE: linear Hardening
   ppr(1) = getrea('UNIAX ',0d0,' Yield Stress .....................')
   ppr(2) = getrea('HARDS ',0d0,' Isotropic Hardening Par. .........')
 CASE (3)  ! IHLUDW: Ludwik Nadai
   ppr(3) = getrea('CONSN ',0d0, '!Isotropic Hardening exponent .....')
   IF(exists('CONS0 ') .AND. exists('EFREF ')) THEN
     ppr(1) = getrea('CONS0 ',1d0,'!Reference Yield Value ............')
     ppr(2) = getrea('EFREF ',0d0,'!EFREF constant  ..................')
   ELSE IF (exists('UNIAX ').AND. &
            exists('EUNIF ').AND.exists('STREN '))THEN
     uniax = getrea('UNIAX ',0d0,'!Elastic limit ....................')
     eunif = getrea('EUNIF ',0d0,'!Limit of Uniform Strain...........')
     strem = getrea('STREN ',0d0,'!Nominal Strength .................')
     ppr(2) = LOG(1d0+eunif)/((strem/uniax)**(1d0/ppr(3))-1d0)
     ppr(1) = uniax/ppr(2)**ppr(3)
   ELSE IF(exists('CONS0 ').AND.exists('UNIAX '))THEN
     ppr(1) = getrea('CONS0 ',1d0,'!Reference Yield Value ............')
     uniax = getrea('UNIAX ',0d0,'!Elastic limit ....................')
     ppr(2) = (uniax/ppr(1))**(1d0/ppr(3))
   END IF
 CASE (4)  ! IHSATU: linear + saturation law
   ppr(1) = getrea('UNIAX ',0d0,'!Yield Stress .....................')
   ppr(2) = getrea('HARDS ',0d0,'!Isotropic Hardening Par. .........')
   ppr(3) = getrea('SATEXP',0d0,'!Saturation Exponent ..............')
   ppr(4) = getrea('REFLST',0d0,'!Residual Flow Stress .............')
 CASE (5)  ! defined by points
   rdcur = .TRUE.  !read curve below
 CASE (6)  ! IHHOLO: Holomon law
   ppr(1) = getrea('CONS0 ',1d0,'!Reference Yield Value ............')
   ppr(3) = getrea('CONSN ',0d0,'!Isotropic Hardening exponent .....')
 CASE (7)  ! IHVOCE: Voce law (saturation)
   ppr(1) = getrea('UNIAX ',0d0,'!Yield Stress .....................')
   ppr(2) = getrea('CONS0 ',1d0,'!Reference Yield Value ............')
   !CONSN or SATEXP is allowed as input data for Voce Saturation Exponent
   IF( exists('CONSN ') )THEN
     ppr(3) = getrea('CONSN ',0d0,'!Saturation Exponent ..............')
   ELSE
     ppr(3) = getrea('SATEXP',0d0,'!Saturation Exponent ..............')
   END IF
   !Check that does not exist both parameters: CONSN and SATEXP. This is to avoid duplicated definition
   IF( exists('CONSN ') .AND. exists('SATEXP') )THEN
     CALL runend ('MATPLA: Both CONSN and SATEXP must not appear')
   END IF
 CASE DEFAULT
   ! nothing
 END SELECT

 ppr(5) = getrea('DEATHT',HUGE(1d0),' End (Deactivation) Time ..........')

 ! second: read kinematic hardening law

 !!If kinematic hardening law not explicitly defined, try to guess it from key-words
 IF( def(5) == 0 )THEN
   IF( exists('KINET '))THEN
     def(5) = 2               !linear kinematic hardening
   ELSE IF( exists('KINMOD').AND.exists('KINSAT'))THEN
     def(5) = 3               !saturation law
   ELSE
     def(5) = 1               !no hardening
   END IF
 END IF

 SELECT CASE (def(5))

 CASE (2)   !KHLINE: linear Hardening
   ppr(6) = getrea('KINET ',0d0,' Kinematic Hardening Par. .........')
 CASE (3)   !KHSATU: saturation law
   ppr(6) = getrea('KINMOD',0d0,' Kinematic Hardening Modulus ......')
   ppr(7) = getrea('KINSAT',0d0,' Kinematic Hardening Saturation ...')
 CASE DEFAULT
   ! nothing
 END SELECT

 !third: read variables associated to plastic anisotropy

 IF( def(3) /= 2 )THEN !not von MISES

   ! check existence of strength ratios

   ! strength ratios
   pars(1) = exists('SXY   ')
   pars(2) = exists('S45   ')
   pars(3) = exists('S90   ')
   pars(4) = exists('BIUNI ')
   npar = COUNT(pars)
   ! initializes
   ! The initialization to 0 implies that no values are read
   ! no default are given, 0 acts as a flag
   ppr(8:11) = 0d0  !s90,biuni,sxy,s45
   IF( npar > 0 )THEN
     IF( pars(3) ) ppr(8)= getrea('S90   ',1d0, &
                                  ' Strength Ratio in local X2 direct.')
     IF( pars(4) ) ppr(9)= getrea('BIUNI ',1d0, &
                                  ' Biaxial strength ratio............')
     IF( pars(1) .AND. pars(2) )THEN
       CALL runend ('INMA02:Both SXY and S45 must not appear')
     ELSE IF( pars(1) )THEN
       ppr(10) = getrea('SXY   ',1d0,' Strength Ratio for local shear    ')
     ELSE IF( pars(2) ) THEN
       ppr(11) = getrea('S45   ',1d0,' Strength Ratio in local skew dir. ')
     END IF

     IF( pars(4) .AND. pars(1) )THEN !BIUNI and SXY
       ppr(11) =2d0/SQRT(3d0/ppr(10)**2+1d0/ppr(8)**2-1d0/ppr(8))
     !ELSE IF( pars(4) .AND. pars(2) ) THEN !BIUNI and S45
     ! ppr(10) =SQRT(3d0)*ppr(9)*ppr(11)/SQRT( 4d0*ppr(9)**2-ppr(11)**2)
     END IF
   END IF

   !plastic anisotropy in terms of Lankford coeffs.
   pars(1) = exists('RM    ')
   pars(2) = exists('R0    ')
   pars(3) = exists('R45   ')
   pars(4) = exists('R90   ')
   ! initializes
   ppr(12:16)  = 0d0

   IF( pars(1) .AND. (pars(2).OR.pars(3).OR.pars(4))) &
      CALL runend ('INMA02:Both RM and R?? must not appear')

   IF( pars(1) )THEN
     ppr(12) = getrea('RM    ',1d0,'!LANKFORD COEFFICIENT (average)....')
   ELSE IF( ANY (pars) )THEN
     ppr(13) = getrea('R0    ',1d0,'!LANKFORD COEFFICIENT IN ROLLING D.')
     ppr(14) = getrea('R45   ',1d0,'!LANKFORD COEFFICIENT DIR. 45 DEG. ')
     ppr(15) = getrea('R90   ',1d0,'!LANKFORD COEFFICIENT DIR. 90 DEG. ')
     !ppr(12) = (ppr(13) + ppr(15) + 2d0*ppr(14))/4d0
   END IF
   IF( exists('MHILL ')) &
     ppr(16) = getrea('MHILL ',2d0,'!HILL exponent ................... ')
   IF (exists('CHILL',i)) THEN
     ppr(16:21) = param(i:i+5)
     ppr(22:27) = param(i:i+5)
   END IF
   IF( ( def(3) == 6 ) .AND. (exists('RB    '))) &
     ppr(16) = getrea('RB    ',1d0,'!LANKFORD COEFFICIENT BIAXIAL  ')     !FR

 END IF

 ! *********** Viscosity parameters
 SELECT CASE (def(6))

 CASE (3)   !ELVIPL: elastic-viscoplastic
   ppr(28) = getrea('VISCO ',0d0,'!Viscosity parameter ...............')
   ppr(29) = getrea('HDEXP ',0d0,'!Strain hardening visco. exponent ..')
   ppr(30) = getrea('RSEXP ',0d0,'!Strain-rate sensivity visco. expo..')
   !QUAD4 viscoplastic formulation same that TR2D (Ponthot-Garino)
   !ppr(16) = ppr(28)
 CASE DEFAULT
   ! nothing
 END SELECT

 IF( rdcur )THEN
   IF(PRESENT(chead) )THEN !read curve for isotropic hardening definition
     CALL listen('MATPLA') !read a line
     IF( exists('CURVED') )THEN     !IF curve defined material read it
       CALL listen('MATPLA') !read line with specification of strain,stress measu
       ! strain measure used in the curve definition
       ! default strain measure: eff.plastic strain
       ! other options:
       ! 'LOGARI' - logarithmic strain (total)
       ! 'ENGINE' - engineering strain
       straim = get_name('STRAIN',defau='EFFPLA',texts='  Strain measure:')
       ! stress measure used in the curve definition
       ! default stress measure: Cauchy stress
       ! other options:
       ! 'NOMINA' - nominal stresses
       stresm = get_name('STRESS',defau='CAUCHY',texts='  Stress measure:')
       ! verification will be added
       CALL rdmpcr(chead,ctail,nc,'STRAIN','STRESS')
     ELSE
       CALL runend ('Curve definition expected for this material')
     END IF
     IF (nc /= 1) CALL runend ('One curve should be given for this material')
     ! conversion of strain and stress measures to default ones (prev. getcrv)
     young = epr(1)
     CALL conv_smeas(chead%val,chead%np,young,straim,stresm)
     ppr(1) = chead%val(2,1) ! Yield Stress
     ppr(2) = chead%val(3,1) ! Hardening
   ELSE
     CALL runend ('IHPOINT not valid for this material')
   END IF
 END IF

 RETURN
 END SUBROUTINE mat_pla
