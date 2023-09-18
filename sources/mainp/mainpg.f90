 PROGRAM Simpact
 !***********************************************************************
 !                        -  S I M P A C T  -
 !                   AN EXPLICIT FINITE ELEMENT PROGRAM
 !
 !          GROUP OF NUMERICAL METHODS IN STRUCTURAL MECHANICS
 !                        CORDOBA, 1994-2013
 !     INTERNATIONAL CENTER FOR NUMERICAL METHODS IN ENGINEERING
 !                        BARCELONA, 1992-2011
 !                           QUANTECH ATZ,
 !                        Barcelona, 1998-2011
 !---------------------------------------------------------------------------
 !     Dr. F.Flores - Main developer and Managment - GMNME- UNC - Argentina
 !---------------------------------------------------------------------------
 !***********************************************************************
!#if UNIX
!   USE ifport, ONLY : BEEPQQ  !call atention in interactive mode
!#else
   USE DFLIB, ONLY : BEEPQQ  !call atention in interactive mode
!#endif
   ! Global data bases
   USE param_db,ONLY : milb, mich                       !Global parameters
   USE ctrl_db, ONLY : ttime ,dtime,endtm,numct,ptype,nstra,text,ndofn,aero   !global control parameters
   USE lispa0,  ONLY : ludat,lures              !
   USE outp_db, ONLY : nreqd,nreqv,nreqa,nreqc,ncdis,toutd,iwrit,time,timed,cpui    !output asked, frequencies
   USE meshmo_db,ONLY: r_meshmo,ntrmstra   !mesh modifications database

   USE export_db, ONLY : export, neset              !export data
   USE gvar_db, ONLY : actchk,static
   USE name_db, ONLY : gen_data_file
   USE introut_sr, ONLY: introut_mng
   IMPLICIT NONE

   !Functions
   CHARACTER(len=mich):: inttoch   !Puts an integer into a string

   !Variables
   CHARACTER(len=milb):: actio='NEW', ptold=''
   INTEGER  (kind=4) :: ircon=0, &  ! counter of restart files
                        nrest,   &  ! restart dumping period
                        kstep
   INTEGER  (kind=4) :: istop=0, &  ! termination/error flag
                        ecdis=0, &  ! control DOF for output
                        dumi!,    &  ! auxiliar
                        !nprocc=1    ! number of processors

   LOGICAL ::           actmsh=.FALSE.   &  ! .TRUE. if mesh modification required in step
                       ,endmsh=.FALSE.   &  ! .TRUE. if mesh modification required at end step
                       ,actrst=.FALSE.      ! .TRUE. if restart file must be written


   INTEGER (kind=4) chnode          ! changes from node Label to internal node number

   INTERFACE
     INCLUDE 'byebye.h'
     INCLUDE 'closef.h'
     INCLUDE 'outdyn.h'
   END INTERFACE

   ! ==============================================================
   !     INITIALIZATION PHASE
   ! ==============================================================

   CALL beginp(actio,nrest,cpui,kstep)

   !  actio
   !     assigned in beginp:
   !     NEW     - new problem
   !     RESTAR  - restart an old problem (no changes or minor)
   !     NSTRA0  - restart an old problem with a new strategy
   !
   !     can be changed later to:
   !     NSTRA1  - a new strategy will involve change of boundary conditions without changing nodes
   !     NSTRA2  - a new strategy will involve change of both boundary conditions and nodes
   !     STOP    - stop the program
   !     MESHMO  - mesh modifications in step

   IF (TRIM(actio) == 'NEW')THEN
     CALL timing(19,1)    !initializes CPU time for the current model
     CALL newprb()        !set output file name root
   ELSE   !OLD problem
     CALL restar(ircon,actmsh,endmsh)     ! read the problem from a re-start file
     ! check if the restart file corresponds to the end of the strategy ==> modify ACTIO
     IF (TRIM(actio) == 'RESTAR' .AND. ttime + dtime/2d0 >= endtm ) actio='NSTRA0'
     IF (TRIM(actio) == 'RESTAR' ) THEN
       CALL wrtpos(1)                       ! write data for post-processing
     END IF
     IF (actmsh) THEN                     ! TRUE if remesh is active before restart
       actio = 'MESHMO'                   ! if restart from remesh in last step before
     END IF
   END IF

   IF (nrest > 0) CALL timerst(.TRUE.,nrest,actrst)   !Initialize time for writing restart file

   DO ! loop over strategies  =====================================

      !           Initialization of a new problem/strategy

      CALL timing(1,1)                    !initializes time for Initial computations
      !IF new strategy with model changes or mesh modifications
      IF (TRIM(actio) /= 'NEW' .AND. TRIM(actio) /= 'RESTAR' ) THEN
        ircon=0     !initializes RESTART counter
        CALL newstr(actio) !Check beginning of a new strategy ==> ACTIO & PTYPE
        CALL timing(18,1)                      !init time for Closing files
        IF(.NOT.actchk)CALL closef(nreqd,nreqv,nreqa,nreqc,numct) !close files for the previous strategy
        CALL timing(18,2)                      !end time for closing files
        IF (TRIM(actio) == 'STOP') EXIT        !EXIT if erroneous data
        IF (endmsh) THEN    ! Mesh modification at end of previous strategy for new one
          CALL timing(11,1)
          WRITE(*,"(/,19x,'*** INITIALIZING NEW STAGE ANALYSIS ***',/)")
          actio = 'NSTRA2'       !nodes are changed
          CALL timing(11,2)
        END IF

      END IF

      SELECT CASE (ptype)

      CASE ("EXPLCT") !EXPLiCiT time integration

        SELECT CASE ( TRIM(actio) )

        CASE ('NEW','NSTRA0','NSTRA2')

          CALL data_inp(actio)

        CASE ('MESHMO')

          ! Mesh modification: 'remeshing' or 'refinement'

          CALL timing(11,1)
          IF (ANY(r_meshmo)) THEN
            CALL remesh(inttoch(ntrmstra,0))
          END IF
          actio = 'NSTRA2'   !nodes are changed
          CALL timing(11,2)

        END SELECT

        !=======  Initial calculations  ===================================

        IF(TRIM(actio) /= 'RESTAR') CALL initial_comp(actio,istop)

        IF( actchk )THEN
          ttime = endtm
          actio = 'NSTRA0'
          CYCLE
        END IF

        !=======  Step-by-step solution  =================================

        CALL timing(1,2)

        IF( ncdis > 10 )THEN    !update control equation for output
          dumi = chnode(ncdis/10)
          ecdis = dumi*10 + MOD(ncdis,10)
        ELSE
          ecdis = ncdis
        END IF

        IF(TRIM(actio) /= 'RESTAR') CALL outdyn(.TRUE., actmsh , ecdis)  !initial values for HISTORY

        IF( static ) THEN !pseudo STATIC analysis

          CALL pstatic( ircon,nrest,istop,actio,ecdis )

        ELSE              !standard DYNAMIC analysis

          CALL dynamic(actio,ircon,nrest,kstep,istop,ecdis,actmsh,endmsh)

        END IF


        IF (actmsh) THEN
          IF (ANY(r_meshmo)) THEN
            WRITE( 6,"(5X,'***  START OF REMESH NUMBER ',A)") TRIM(inttoch(ntrmstra+1,0))
            WRITE(55,"(5X,'***  REMESH NUMBER ',A)",ERR=9999) TRIM(inttoch(ntrmstra+1,0))
          END IF
        ELSE
          WRITE( 6,"(5X,'***  END OF STAGE ',A,': ',A,/)") TRIM(inttoch(nstra,0)),text
          WRITE(55,"(5X,'***  END OF STAGE ',A,': ',A,/)",ERR=9999) TRIM(inttoch(nstra,0)),text

          IF( neset > 0 ) CALL export( ) !export selected element sets

          IF (endmsh) THEN ! mesh modification at end strategy
          END IF
        END IF

      CASE ("TRAROT")  !  *******  POSITIONING  *******
        CALL trarot(ecdis,ircon)

      END SELECT

      IF (istop /= 0)EXIT

      ptold = ptype
      IF (actmsh) THEN
         actio = 'MESHMO'
      ELSE
         actio = 'NSTRA0'
      END IF

   END DO  ! end loop over strategies

!=======================================================================
!  ******   ENDING EXECUTION
!=======================================================================
   IF(istop /= 0) THEN
      WRITE(* ,"(/,20x,'*** END PROGRAM - PROBLEMS DETECTED ***',/)")
      WRITE(55,"(/,20x,'*** END PROGRAM - PROBLEMS DETECTED ***',/)",ERR=9999)
      WRITE(55,"(/,22x,'*** ERRONEOUS RESULTS POSSIBLY ***',/)",ERR=9999)
   ELSE
      WRITE(* ,"(/,23x,'*** NORMAL END OF EXECUTION ***',/)")
      WRITE(55,"(/,23x,'*** NORMAL END OF EXECUTION ***',/)",ERR=9999)
   END IF
   CALL timing(19,2)
   CALL byebye(time)
   !CALL debtime(timed(1))
   IF( gen_data_file )THEN
     CLOSE(ludat,STATUS='delete')  !PROGRAM.dat file
   ELSE
     CLOSE(ludat,STATUS='keep')    !Users data file
   END IF
   CLOSE(lures,STATUS='KEEP') !
   CLOSE (55)                 !extended report file
   CLOSE (58)                 !debug file for developers
   
   IF( aero )THEN
     CALL introut_mng('end')
   ENDIF

   CALL BEEPQQ(5000,50)  !call atention
   STOP
 ! -----------------------
 9999 CALL runen2('ERROR WHILE WRITING "REP" FILE.')

 END PROGRAM Simpact
