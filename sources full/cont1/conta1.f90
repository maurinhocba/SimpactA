 SUBROUTINE conta1(itask,dtcal,ttime,iwrit,velnp,maxve)

  !     main  routine for tube-skin 3D contact

 USE cont1_db
 IMPLICIT NONE

   !        Dummy arguments
                                 !task to perform
   CHARACTER(len=*),INTENT(IN):: itask
   INTEGER(kind=4),INTENT(IN):: iwrit,maxve
   REAL(kind=8),INTENT(IN),OPTIONAL:: dtcal,ttime,velnp(:,:)
   !        Local variables
   TYPE (pair1_db), POINTER :: pair
   INTEGER (kind=4) :: ipair,i
   REAL (kind=8) ::  dtime, auxil
   REAL (kind=8), SAVE :: timpr = 0d0

   !....  PERFORM SELECT TASK

   SELECT CASE (TRIM(itask))

   CASE ('NEW','NSTRA0','NSTRA1','NSTRA2')   !Read Input Data
       !from npo_db INTENT(IN) :: npoin,coord,label
       ! .. Initializes data base and input element DATA
       CALL cinpu1(maxve,iwrit)

   CASE ('LUMASS')
     !  compute mass
     CALL surms1(emass)

   CASE ('FORCES')  !Performs contact search and computes contact forces
     !from npo_db INTENT(IN) :: npoin,coora,label,emass
     !            INTENT(IN OUT) :: fcont
     IF( ctime < 0d0 )ctime = dtcal      !first step only
     dtime = ctime                       !contact dtime from database
     IF(dtime == 0d0)dtime = dtcal       !computation time

     ! ....  Perform Contact Searching & Compute contact forces
     CALL celmn1(ttime,dtime,coora,emass,fcont,velnp)
     timpr = timpr + dtime               !increase elapsed time
   CASE ('DUMPIN')   !write data to a re-start file
     CALL cdump1()

   CASE ('RESTAR')   !read data from a re-start file
     CALL crest1()
     ALLOCATE ( surtf(3,npair) )              !get memory for total forces
     surtf = 0d0                              !initializes

   CASE ('UPDLON')   !Modifies internal numeration
     CALL cupdl1( )

   CASE ('OUTDY1')   !Writes contact forces between pairs
     dtime = ctime                       !contact dtime from database
     IF(dtime == 0d0)dtime = dtcal       !computation time
     WRITE(41,ERR=9999) ttime            !control variable
     pair => headp
     DO ipair=1,npair                      !for each pair
       IF(timpr > 0d0 )THEN
         WRITE(41,ERR=9999) surtf(1:3,ipair)/timpr     !average contact forces
       ELSE
         WRITE(41,ERR=9999) surtf(1:3,ipair)           !average contact forces
       END IF
       IF( pair%press ) pair%presn = pair%presn/(timpr+dtime)*dtime
       pair => pair%next
     END DO
     surtf = 0d0                             !initializes for next period
     timpr = 0d0                             !initializes elapsed time
   CASE ('OUTDY2')   !Writes contact forces and gaps for nodes
     dtime = ctime                       !contact dtime from database
     IF(dtime == 0d0)dtime = dtcal       !computation time
     auxil = timpr+dtime
     IF( auxil == 0d0 ) auxil = 1d0
     pair => headp
     DO ipair=1,npair                      !for each pair
       IF( pair%press ) WRITE(44,ERR=9999) (pair%presn(i)/auxil,i=1,pair%ncnod)
       pair => pair%next
     END DO
   CASE ('INITIA')   !Control Initial penetrations
     ! ....  Perform Contact Searching & Output a report
     CALL celmn1p(ttime,coora)
   CASE ('WRTSUR')   !Writes contact surfaces for Post-processing
     ! OPEN files for postprocess, detect surface to use
     CALL prsur1(iwrit)
   END SELECT

 RETURN
 9999 CALL runen2('')

 END SUBROUTINE conta1
