 SUBROUTINE dumpin(ircon,endst,actmsh,endmsh)

! writes restart file

 USE param_db,ONLY: mlin
 USE c_input,ONLY: openfi
 USE ctrl_db
 USE curv_db
 USE damp_db
 USE esets_db,ONLY: eset, msets, esets, nelms, dump_esets
 USE ifx_db
 USE kinc_db, ONLY: dump_kinc, naris, nvelr
 USE mat_dba
 USE nsets_db, ONLY: dump_nsdb
 USE npo_db
 USE outp_db, ONLY: dump_outp
 USE meshmo_db
 USE gvar_db,ONLY: nrst,static
 USE ift_db,ONLY: dump_ift
 USE flc_db,ONLY: dump_flc
 USE cms_db,ONLY: dump_cm
 USE loa_db,ONLY: dump_lo
 USE export_db,ONLY: dump_exp_data
 USE sms_db,ONLY: dump_sms_data
 USE static_db,ONLY: dump_static
 USE nsld_db,ONLY: nnsld

 IMPLICIT NONE

   ! Dummy variables
   INTEGER(kind=4),INTENT(IN OUT):: ircon
   LOGICAL,INTENT(IN):: endst,  & !.TRUE. if the end of the strategy has been reached
                        actmsh, & !.TRUE. if last strategy modified the mesh
                        endmsh    !.TRUE. if end estrategy mesh modification is required
   ! Local variables
   CHARACTER(len=mlin) :: namei
   CHARACTER(len=1),SAVE ::  ni(0:9)=(/'0','1','2','3','4','5','6','7','8','9'/)
   INTEGER (kind=4) :: ircod,j2,j3,j4,flag
   REAL (kind=8) :: cpuac

   INTERFACE
     INCLUDE 'contac.h'
     INCLUDE 'elemnt.h'
   END INTERFACE

   !       ircon defines the file where to WRITE
   ircon = ircon+1
   !       ircod defines the old file to delete (if exist)
   ircod = ircon-nrst  !keeps latest 'nrst' files for restart

   j2 = MOD(ircon,1000)/100
   j3 = MOD(ircon,100)/10
   j4 = MOD(ircon,10)
   namei = '.'//ni(j2)//ni(j3)//ni(j4)   !file name
   IF (endst) THEN  !file corresponds to the end of the strategy
     flag = 0  ! file corresponds to end of the strategy
   ELSE
     flag = 1  ! file corresponds to period NREST
   END IF

   CALL timuse(cpuac)    !compute actual clock time
   CALL openfi(nfile=50,rest=namei,flag=flag,oflag=flag)  !open file to write in

   WRITE(50,ERR=9999) ircon, cpuac
   WRITE(50,ERR=9999) actmsh, endmsh

   CALL dump_ctrl        !dumps CTRL_DB
   CALL dump_esets       !dumps ESETS_DB
   CALL dump_nsdb  ! dumps nodal sets database
   CALL dump_kinc(nnsld)  ! dumps kinematic constrains parameters
   CALL dump_npo   ! dumps nodal data base
   CALL dump_outp  ! dumps output control parameters

   CALL dump_curv    ! dumps curve database
   CALL dump_flc()   ! dumps FLC curve database
   CALL dump_mate    ! dumps material database

   ! dump element databases
   CALL elemnt ('DUMPIN')

   WRITE(50,ERR=9999) 'ENDSET'

   CALL dump_kc (naris, nvelr, ndofn)
   CALL dump_ifx
   CALL dump_cm (ndofn)
   CALL dump_sms_data ()
   CALL dump_lo (ndime,ndofn,nload)
   CALL dump_damp (neq)
   CALL dump_meshmo

   IF(ANY(cactive)) CALL contac('DUMPIN',0)

   IF(therm) THEN
     CALL dump_ift
   END IF
   CALL dump_exp_data
   WRITE(50,ERR=9999) static
   IF( static )CALL dump_static (neq,ndime,npoin)

   CLOSE(50,STATUS='keep')
   CALL flushf(50)

   IF (flag /= -1) THEN
     CALL manage_last_file_names(1)   !Stores the name of last restart file
     IF (ircod > 0 .AND. flag /= -1) THEN
       j2 = MOD(ircod,1000)/100
       j3 = MOD(ircod,100)/10
       j4 = MOD(ircod,10)
       namei = '.'//ni(j2)//ni(j3)//ni(j4)
       CALL del_old_files(flag=3,rest=namei)
     END IF
   END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dumpin
