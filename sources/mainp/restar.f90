 SUBROUTINE restar (ircon,actmsh,endmsh)

!  read a restart file

 USE ctrl_db
 USE curv_db
 USE damp_db
 USE esets_db,ONLY: eset, msets, nelms, esets, rest_esets
 USE ifx_db
 USE kinc_db, ONLY: rest_kinc, naris, nvelr
 USE c_input
 USE mat_dba
 USE nsets_db, ONLY: rest_nsdb
 USE npo_db
 USE outp_db, ONLY : rest_outp,time,cpui
 USE gvar_db,ONLY: static
 USE meshmo_db
 USE param_db,ONLY: mlen
 USE name_db,ONLY: output, namr, rsfprj, rsfstg
 USE ift_db,ONLY: rest_ift
 USE flc_db,ONLY: rest_flc
 USE cms_db,ONLY: rest_cm
 USE loa_db,ONLY: rest_lo
 USE export_db,ONLY: rest_exp_data
 USE sms_db,ONLY: rest_sms_data
 USE static_db,ONLY: rest_static
 USE nsld_db,ONLY: nnsld
 USE structuresa,  ONLY: restar_structuresa
 USE arraysa,      ONLY: restar_arraysa
 USE aeroout_db,   ONLY: restar_aeroout
 USE wind_db,      ONLY: restar_win
!USE gener_db,     ONLY: restar_gen
 USE inter_db,     ONLY: restar_inter
 USE introut_db,   ONLY: restar_introut
!USE mastrnods_db, ONLY: restar_mastrnods
 USE struout_db,   ONLY: restar_struout
 USE inter_sr,     ONLY: introut_mng

 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4),INTENT(OUT) :: ircon
   LOGICAL,INTENT(OUT):: actmsh, endmsh
   !--- Local variables
   INTEGER(kind=4) :: k,kstra,lrn,ist
   REAL(kind=8) :: cpuol, cpuac

   INTERFACE
     INCLUDE 'contac.h'
   END INTERFACE

   cpuol = cpui !stores initial CPU time of present run

   CALL openfi(nfile=51,rest=namr)   !open restart file (NAMR)

   READ(51) ircon, cpuac
   READ(51) actmsh, endmsh   !Flags indicating if a mesh change was performed as the last task

   CALL rest_ctrl    !Read control parameters
   CALL rest_esets   !Read esets db
   CALL rest_nsdb    !restores nodal sets database
   CALL rest_kinc(nnsld)    !restores kinematic constrains parameters
   CALL rest_npo     !restores NPO arrays
   CALL rest_outp    !reads output control parameters


   !Total Process Time
   time(19) = time(19) + cpuac - time(39)
   time(39) = cpuol
   cpui = cpuol - time(19)

   CALL rest_curv    ! restores curve database
   CALL rest_flc()   ! restores FLC curve database
   CALL rest_mate
   CALL rest_el

   CALL rest_kc (naris, nvelr, ndofn, ndime)
   CALL rest_ifx
   CALL rest_cm (nrotd)
   CALL rest_sms_data (ndime)
   CALL rest_lo (ndime,nrotd,nload)
   CALL rest_damp (neq)
   CALL rest_meshmo()

   lrn = lfile(namr)                       !Lenght of the original name project
   IF (TRIM(output)/=namr(1:lrn)) ircon=0  !Restart counter reset
   CALL defnam(output,nstra)   !Determina el nombre del fichero de salida
   lrn = lfile(output)
   rsfprj = output(1:lrn)                   !Name of the project  (eg. data file: <rsfprj>.dat)
   rsfstg = output(lrn+1:LEN_TRIM(output))  !Name of the strategy (eg. output file: <rsfprj><rsfstg>.f16)

   !IF(numct > 0)   CALL contac('RESTAR')
   IF(ANY(cactive))   CALL contac('RESTAR',0)
   !      Top and bottom surfaces used by contact and coupled problems in shells
   IF( bottom ) ALLOCATE (coorb(ndime,npoin))      !get memory for contact bottom surface
   IF( top    ) ALLOCATE (coort(ndime,npoin))      !get memory for contact top surface
   IF( top .OR. bottom ) ALLOCATE ( ifact(npoin) ) !get memory for contact auxiliar array

   IF(therm) THEN
     CALL rest_ift
   END IF
   CALL rest_exp_data
   READ (51) static
   IF( static )CALL rest_static (neq,ndime,npoin)
   
   IF(aero)THEN
     CALL restar_structuresa
     CALL restar_arraysa
     CALL restar_aeroout
     CALL restar_win
     !CALL restar_gen
     CALL restar_inter
     CALL restar_introut
     !CALL restar_mastrnods
     CALL restar_struout
     CALL introut_mng('openFiles','RESTAR')
   ENDIF

   CLOSE(51,STATUS='keep')

   ! here the search of the end of the last strategy
   ! after this change there must be only one data file

   kstra = 0
   DO
     READ(ludat,"(A)") card
     k = LEN_TRIM(card)
     IF( k > 0 )THEN
       CALL upcase(card(1:k))
       IF( INDEX(card(1:k),'STR') /= 0)THEN
         BACKSPACE ludat  ! do not replace by   backs = .TRUE.
         CALL listen ('RESTAR')
         IF (exists('ENDSTR')) kstra = kstra + 1
         IF (kstra == nstra) EXIT
       END IF
     END IF
   ENDDO

   CALL manage_last_file_names(2)

   WRITE(*,"('+ PROGRAM SUCCESFULLY RESTARTED')")
   WRITE(55,"('+ PROGRAM SUCCESFULLY RESTARTED')",IOSTAT=ist)
   IF (ist /= 0) CALL runen2('ERROR WHILE WRITING "DEB" FILE.')

 RETURN
 END SUBROUTINE restar
