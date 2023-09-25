 SUBROUTINE contol(actio)
 !******************************************************************
 !*** READ control DATA
 !******************************************************************
 USE param_db,ONLY: mnam, mttl
 USE ctrl_db,ONLY: ndime, ndofn, nrotd, neulr, ntype, nacce, ifunc, ifixd, dtrec, endtm,  &
                   begtm, ncdlt, mscal, dtscal, dtuser, text, dtime,ttime, vefac,  &
                   lumped, itemp, therm, ndoft, aero,            &
                   tdtuse, tdtime, tdtsca, tmscal, tref, tscal, ascal
 USE c_input,ONLY: listen, exists, ludat, lures, get_name, getrea, getint,  &
                   card, backs, words
 USE outp_db, ONLY: energ, postype, toutp, toutd, ncdis, iener, cname, iwrit
 USE meshmo_db,ONLY: r_nrfset
 USE npo_db, ONLY : acceg
 USE esets_db, ONLY : nrenu
 USE gvar_db, ONLY : actchk,static
 USE static_db, ONLY : nides,miter,ittol,tol,sstep,pcont,nperi,lovec,nss,nds, &
                       sstop,autod,pmin,pmax,ipred,tolm,fdamp,spback
 IMPLICIT NONE

  !--- Dummy variables
  CHARACTER(len=*),INTENT(IN):: actio
  !--- Local variables
  INTEGER(kind=4):: i, j, point, dof
  LOGICAL done(5), dout,tout,cout

  INTERFACE
    SUBROUTINE readtp(toutp,t_d,task)
      IMPLICIT NONE
      CHARACTER(len=1), INTENT(IN) :: t_d
      INTEGER (kind=4), INTENT(IN) :: task
      REAL (kind=8), POINTER :: toutp(:)
    END SUBROUTINE readtp
  END INTERFACE

  !*** READ the first DATA card, and echo it immediately.
  IF (TRIM(actio) == 'NEW') THEN
    done = .FALSE.      !initializes all tasks to FALSE.
    !  assign default values for some parameters
    dtuser = 0d0        !time increment
    ncdlt = 0           !do not recompute
    dtscal = 0.75d0     !time increment scale factor
    mscal  = 1d0        !mass scale factor
    iwrit  = 1          !flag to echo data and results writing
    iener  = 0          !flat to compute energy parameters
    !        default values for heat transfer
    tscal  = 1.0d0      !time scale factor for thermal analysis
    therm  = .FALSE.    !thermo-mechanical analysis
    tdtuse = 0d0        !time increment for heat transfer
    tdtsca = 0.75d0     !time increment scale factor for heat transfer
    tmscal = 1d0        !mass scale factor for heat transfer
    ascal = .FALSE.     !don't auto-compute mass scale factor for thermo-mech coupling
    tref   = 293d0      !reference temperature (20 Celsius)
    itemp  = .FALSE.    !there are temperatures in the model
    ndoft  = 0          !no temperature DOFs
  ELSE
    done = .TRUE.       !initializes all flags to TRUE
    cout = .FALSE.      !only OUTPUT control is compulsory
  END IF

  CALL listen('CONTOL')            !read a line
  DO
    IF (exists('ENDCON')) EXIT     !end line found, exit loop
    ! C O N T R O L    D A T A   (new problem only)
    IF (exists('CONTRO')) THEN     !Set of control variables
      IF (TRIM(actio) /= 'NEW') CALL runend('CONTOL: CONTROL Params cannot chang')  !check
      WRITE (lures,"(/,'  C O N T R O L   P A R A M E T E R S ',/)",ERR=9999)
      ndime=getint('NDIME ',3,' Number of coordinate Dimension ...')
      IF( ndime == 2 )ntype=getint('NTYPE ',3,' Problem type for 2 Dimensional ...')
      neulr=getint('NEULR ',0,' Nodal system code 0:NO 1:YES .....')
      ndofn=getint('NDOFN ',ndime,' Number of Degrees Of Freedom .....')
      nrotd = MIN(ndofn,6)
      IF( exists('CONSIS') )THEN
        lumped = .FALSE.  !consistent mass matrix
        WRITE (lures, "(/,'      Consistent Mass Matrix will be used',/)",ERR=9999)
        nrenu=getint('NRENU ',0,' Node Renumbering flag ............')
      END IF
      WRITE (lures, "(/,'Type of analysis:')",ERR=9999)
      IF (exists('COUPLE')) THEN
        WRITE (lures, "(/,'      coupled thermo-mechanical',/)",ERR=9999)
        therm = .TRUE.
      ELSE
        therm = .FALSE.
        WRITE (lures, "(/,'      mechanical only',/)",ERR=9999)
      END IF
      IF (exists('AERO')) THEN
        WRITE (lures, "(/,'      Aero-elastical problem   ',/)",ERR=9999)
        aero  = .TRUE.
      ELSE
        aero  = .FALSE.
      END IF
      ndoft=getint('NDOFT ',1,' Number of Temperature DOFs..........')
      done(1)=.TRUE.    ! set first task as done
    END IF
    !  **** S T E P   T I M E   D A T A
    IF (exists('TIME  ') .OR. exists('STATIC')) THEN     !time control parameters
      WRITE (lures,"(/,'  T I M E   P A R A M E T E R S ',/)",ERR=9999)
      IF (TRIM(actio) == 'NEW') THEN     !compulsory variable for NEW problems
        begtm = 0d0 ! initial time of the new strategy
        endtm =getrea('ENDTM ',0d0,'!Final Time .......................')
      ELSE                         !optional for regular strategies
        begtm = ttime ! initial time of the new strategy (=final of the last
        endtm =getrea('ENDTM ',endtm,' Final Time .......................')
      END IF
      dtscal=getrea('DTSCAL',dtscal,' SCAle fact. for DT in auto mode...')
      IF (dtscal == 0d0) dtscal=0.75d0
      IF (dtscal > 0.8d0) THEN
        WRITE(*,"(/,' Warning! Scaling factor DTSCAL might be too HIGH')")
        WRITE(55,"(/,' WARNING! Scaling factor DTSCAL might be too HIGH',/)",ERR=9999)
      END IF
      ncdlt = getint('NCDLT ',ncdlt,' No of steps between DT computat..')
      IF (exists('STATIC')) THEN     !static control parameters
        static = .TRUE.
        spback= exists('SPRING')      !elastic springback type
        miter =getint('MITER ',10000,' Maximum number of iterations .....')
        ittol =getint('ITTOL ',1    ,' Tolerance type ...................')
        tol   =getrea('TOLER ',0.01d0,' Convergence Tolerance ............')
        tolm  =getrea('TOLMAX',5d0*tol,' Maximum Convergence Tolerance ...')
        IF( spback ) THEN
          sstep = 1
          fdamp = 1d0
          sstop = .TRUE.
          autod = .FALSE.
          nides = miter/100
          nperi = miter/10
          ipred = 1 ; lovec = 1 ; nss = 1 ; nds = 1
        ELSE
          sstep =getint('SSTEP ',10   ,' Number of steps for the stage ....')
          fdamp =getrea('FACTOR',1.4d0,' Damping factor between steps......')
          sstop = exists('STOP')        !stop if no convergence
          autod = exists('AUTODA')      !use automatic damping
          nides =getint('NIDES ',100  ,' Number of steps to apply displ....')
          nperi =getint('NPERI ',miter/10,' Minimum Number of iteration   ....')
          ipred =getint('IPRED',1,     ' Prediction order .................')
          lovec =getint('LVSET ',0    ,' Load (+) or Velocity (-) control .')
          nss   =getint('NSS   ',2    ,' Number of subdivisions of veloc. .')
          nds   =getint('NDS   ',5    ,' Number of displacements steps ....')
          pmin  =getint('PMIN  ',   1 ,' Minimum Damping period ...........')
          pmax  =getint('PMAX  ',miter,' Maximum Damping period ...........')
        END IF
      ELSE
        dtuser=getrea('DTUSER',dtuser,' Time Increment = 0.-> Auto calc...')
        IF( dtuser == 0d0 .AND. ncdlt == 0 ) ncdlt = 100    !frequency of dtime re-evaluation
        mscal =getrea('MSCALE',mscal,' Mass SCALing factor  .............')
        tscal =getrea('TSCALE',1d0,' Time SCALing factor for therm.ana.')
        IF (exists('DTUSER') .AND. dtuser > 0d0 )THEN
          IF( dtime < dtuser .AND. TRIM(actio) /= 'NEW')THEN
             WRITE (lures,"(/,'  W A R N I N G   DTIME modified by user',/ &
                        &     '  DTUSER:',e12.4,' Old_DTIME:',e12.4)",ERR=9999)dtuser,dtime
             WRITE (*,"(/,'  W A R N I N G   DTIME modified by user',/  &
                        &     '  DTUSER:',e12.4,' Old_DTIME:',e12.4)",ERR=9999)dtuser,dtime
          END IF
          dtime = dtuser
        END IF
        vefac=getrea('VEFAC',vefac,' Velocity smoothing factor ........')
        static = .FALSE.
      END IF
      done(2) = .TRUE.  ! set second task as done
    END IF
    !  **** O U T P U T   C O N T R O L
    IF (exists('OUTPUT')) THEN   ! control parameters for output results
      WRITE (lures,"(/,'  O U T P U T   P A R A M E T E R S ',/)",ERR=9999)
      point = ncdis/10           ! node point
      dof   = MOD(ncdis,10)      ! DOF
      IF (dof == 0) dof=ndime   ! arbitrary set to Y (2D) or Z (3D)
      point = getint('POINT ',point,' Node for Output Control ..........')
      dof   = getint('DOF   ',dof  ,' Dof  for Output Control ..........')
      ncdis = point*10+dof       ! keep in one variable
      ! see if new output control exists
      j = 0   !initializes number of pos-types
      IF (exists('COUTD ') .OR. exists('COUTP ')) THEN   !curve controlled
        cout = .TRUE.                                    !set to true
        j = j + 1                                        !increase counter
        IF (exists('CURVE ',i )) THEN                    !name of the curve exists
          cname = get_name(posin=i,stype='CURV')         !get curve name
        ELSE
          CALL runend(' CURVE key word is compulsory')   !unnamed curve
        END IF
      ELSE
        cout = .FALSE.                                   !no curve data
      END IF
      IF (exists('TOUTD ') .OR. exists('TOUTP ')) THEN   !time controlled
        tout = .TRUE.                                    !set to true
        j = j + 1                                        !increase counter
      ELSE
        tout = .FALSE.                                   !no time  data
      END IF
      IF (exists('DOUTD ') .OR. exists ('DOUTP ')) THEN  !displacement controlled
        dout = .TRUE.                                    !set to true
        j = j + 1                                        !increase counter
        IF (ncdis < 10) CALL runend('CONTOL: check POINT and DOF')
      ELSE
        dout = .FALSE.                                   !no displacement data
      END IF

      IF (j > 1) CALL runend('CONTOL: check output control pars ')

      IF (TRIM(actio) == 'NEW') THEN  !for a new problem
        IF (cout) THEN             ! if a curve data will be used as control
          ncdis = -1               !use negative as a flag
          toutd = getrea('COUTD ',0d0,'!History Output function Interval ..')  !frequency measured as curve incr.
          postype = 'C'                   !set to pos-type CURVE
          CALL readtp(toutp,postype,1)    !read frequency for FULL output
        ELSE IF (dout) THEN      ! if a node displ. will be used as control
          toutd = getrea('DOUTD ',0d0,'!History Output Travel Interval ...')   ! frequency measured as displ. incr.
          postype = 'D'                   ! set to pos-type DISPLACEMENT
          CALL readtp(toutp,postype,1)
        ELSE                     ! if TIME will be used as control
          toutd =getrea('TOUTD ',0d0,'!History Output Interval  .........')  ! time frequency for selected output
          postype = 'T'                   ! set to pos-type TIME
          CALL readtp(toutp,postype,1)    ! read frequency for FULL output
          ncdis = 0
        END IF
      ELSE                       ! for an existing problem
        IF (cout) THEN           ! if a curve data will be used as control
          ncdis = -1
          toutd = getrea('COUTD ',toutd,' History Output Function Interval .')  !frequency measured as curve incr.
          postype = 'C'                   !set to pos-type CURVE
          CALL readtp(toutp,postype,2)    ! read frequency for FULL output
        ELSE IF (dout) THEN      ! if a node displ. will be used as control
          toutd = getrea('DOUTD ',toutd,' History Output Travel Interval ...')  ! frequency measured as displ. incr.
          postype = 'D'                   ! set to pos-type DISPLACEMENT
          CALL readtp(toutp,postype,2)    ! read frequency for FULL output
        ELSE                     ! if TIME will be used as control
          ncdis = 0
          toutd = getrea('TOUTD ',toutd,' History Output Interval  .........')  ! time frequency for selected output
          postype = 'T'                   ! set to pos-type TIME
          CALL readtp(toutp,postype,2)    ! read frequency for FULL output
        END IF
      END IF

      IF( toutd == 0d0 )toutd = 1d0
      iwrit = getint('IWRIT ',iwrit,' Print-out code -      0:NO 1:Yes .')   !write flag
      iener = getint('IENER ',iener,' Calculation of Energy 0:NO 1:Yes .')   !energy computation flag
      done(3) = .TRUE.  ! set third task as done
      IF( static )THEN
        point = getint('CPOINT',point,' Node for Convergence Output.......')
        dof   = getint('CDOF  ',dof  ,' Dof  for Convergence Output.......')
        pcont = point*10+dof       ! keep in one variable
      END IF
      ! curve data position is set in NCDIS, this value is assigned in EQNUMS
      IF( cout ) CALL rdcurv('COUTP', cname )
    END IF
    ! G R O U N D   A C C E L E R A T I O N   D A T A
    IF (exists('ACCELE')) THEN
      !STOP 'ACCELEROGRAMs not active presently'
      WRITE (lures,"(/,'  A C C E L E R O G R A M  PARAMETERS ',/)",ERR=9999)
      dtrec = getrea('DTREC ',0d0,'!Time increment of Accelerogram....')
      nacce = getint('NACCE ',1,'!Points in the Accelerogram .......')
      ifixd = getint('IFIXD ',0,' Direction of excitation ......... ')
      ifunc = getint('IFUNC ',1,' Activate  0:Yes  /=0:No ......... ')
      done(4) = .TRUE.  ! set fourth task as done
    END IF
    ! C O U P L E D   T H E R M O - M E C H   D A T A
    IF (exists('THERM ')) THEN     !heat transfer control parameters
      WRITE (lures,"(/,'  H E A T   T R A N S F E R   P A R A M E T E R S',/)",ERR=9999)
      tdtuse=getrea('TDTUSE',tdtuse,' Time Increment = 0.-> Auto calc...')
      tdtsca=getrea('TDTSCA',tdtsca,' SCAle fact. for DT in auto mode...')
      tmscal=getrea('TMSCAL',tmscal,' Capacity SCALing factor ..........')
      IF(exists('AUTOSC')) &
        ascal = .TRUE. !automatic evaluation of MASS MATRIX SCALING FACTOR
      tref  =getrea('TREF  ',tref,' Reference temperature ............')
      IF(tdtsca == 0.) tdtsca = 0.75d0
      IF(tdtsca > 0.8d0)THEN
        WRITE(* ,"(/,' Warning! Scaling factor TDTSCA might be too HIGH')")
        WRITE(55,"(/,' WARNING! Scaling factor TDTSCA might be too HIGH',/)",ERR=9999)
      END IF
      ! curve data position is set in TCURV, this value is assigned in TEQNUMS
      itemp = .TRUE.
      done(5) = .TRUE.  ! set five task as done
    END IF

    CALL listen('CONTOL')
    IF (ALL(done(1:5)) .AND. TRIM(actio) == 'NEW') THEN
      IF (.NOT.exists('ENDCON')) backs = .TRUE.
      EXIT    !exit if all tasks done and no END_CONTROL card exists
    END IF
  END DO
  !-------       Check all information exists
  IF (.NOT.done(1)) CALL runend('CONTOL:  Specify CONTROL Parameters')
  IF (.NOT.done(2)) CALL runend('CONTOL:  Specify TIME Parameters   ')
  IF (.NOT.done(3)) CALL runend('CONTOL:  Specify OUTPUT Parameters ')
  IF (.NOT.done(4)) THEN  !ground accelerations not read
    dtrec = 0d0         !
    nacce = 1
    ifixd = 0
    ifunc = 1
  END IF
  IF (therm .AND. .NOT.done(5) )THEN
    CALL runend('CONTOL:  Specify THERMAL Parameters')
  ELSE IF (.NOT.therm )THEN
    WRITE (lures, "(/,'Type of analysis: mechanical',/)",ERR=9999)
  END IF
  !-- deal with arrays ACCEG and ENERG
  ALLOCATE(acceg(nacce,ndime), energ(6*iener))
  acceg = 0d0  ! to avoid problems in the debugger version
  energ = 0d0  ! to avoid problems in the debugger version
  !-------
  CALL listen('CONTOL')     !refinament control data
  IF (TRIM(actio) == 'NEW') THEN      !for a new problem
    r_nrfset = 0   !initializes remesh modifications flag
  END IF
  ! M E S H   M O D I F I C A T I O N S   D A T A
  IF (exists('MESHMO')) THEN   ! mesh modifications data exist
    CALL rdmeshmo()
  ELSE                  ! no data
    backs = .TRUE.
  END IF

  IF (ndime == 3 .AND. neulr > 0) neulr=9    ! Local system in 3D defined by a matrix
  IF(actchk) iwrit = 0

RETURN
 9999 CALL runen2('')
END SUBROUTINE contol

SUBROUTINE readtp(toutp,t_d,task)
!-----------------------------------------------------------------------
! read control parameters for full output
!-----------------------------------------------------------------------
USE param_db,ONLY: midn, mlin
USE c_input,ONLY: exists, words, lures, param, getrea
IMPLICIT NONE

  !--- Dummy variables
  CHARACTER(len=1),INTENT(IN):: t_d    ! Time or Displacement control
  INTEGER(kind=4),INTENT(IN):: task   ! 1=NEW  2=OLD
  REAL(kind=8),POINTER:: toutp(:)     ! INTENT(IN OUT)

  !--- Local variables
  CHARACTER(len=midn):: name
  CHARACTER(len=mlin):: label
  INTEGER (kind=4) :: i,j,k
  REAL (kind=8) :: defau, first

  name = t_d // 'OUTP '           !TOUTP or DOUTP to read from data file
  IF (t_d == 'C') THEN            ! curve controlled output
    label = ' Complete Output Curve IntervaL ...'
  ELSE IF (t_d == 'T') THEN       ! time controlled output
    label = ' Complete Output IntervaL .........'
  ELSE                            ! displacement controlled output
    label = ' Complete Output Travel IntervaL ..'
  END IF

  IF (task == 1) THEN            ! for a NEW problem
    label(1:1) = '!'             ! parameters are COMPULSORY
    defau = 1d0                  ! no value as default
  ELSE                           ! for old problems
    defau = toutp(2)             ! use previous value as default
  END IF

  IF (exists(TRIM(name),i)) THEN  ! if ?OUTP exists
    j = i
    DO
      IF ((LEN_TRIM(words(j+1)) > 0) .OR. (param(j+1) == 0d0)) EXIT
      j = j + 1
    END DO
    first = getrea(name,defau,label)
    k = j - i + 2                      !array dimension
    IF (task == 2) DEALLOCATE(toutp)   !for an OLD problem release space
    ALLOCATE(toutp(k))                 !get space for data
    toutp(1) = REAL(k-1,8)             !number of prescribed values
    toutp(2:k) = param(i:j)            !prescribed values for output
    IF (k > 2) WRITE(lures,"(5e15.5)",ERR=9999) toutp(3:k)
  END IF
  IF( toutp(2) == 0d0 ) toutp(2) = 1d0
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE readtp
