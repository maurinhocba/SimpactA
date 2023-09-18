 SUBROUTINE conta3(itask,dtcal,ttime,iwrit,velnp,maxve)

 !     main contac routine (ALGORITHM 3)

 USE ctrl_db, ONLY:  bottom,top ,tdtime
 USE npo_db, ONLY : coord,emass,fcont,coorb,coort
 USE cont3_db
 USE gvar_db, ONLY : static,updiv

 IMPLICIT NONE

 !        variable for statistical record
 INTEGER (kind=4), SAVE :: nn(41,3)
 REAL (kind=8) :: td
 !        Dummy arguments
                                !task to perform
 CHARACTER(len=*),INTENT(IN):: itask
 INTEGER(kind=4),INTENT(IN) :: iwrit
 INTEGER(kind=4),INTENT(IN), OPTIONAL :: maxve
 REAL(kind=8),INTENT(IN),OPTIONAL ::dtcal,ttime,velnp(:,:)
 !        Local variables
 TYPE (pair3_db), POINTER :: pair
 INTEGER (kind=4) :: ipair,i
 REAL (kind=8) :: disma, dtime, auxil
 REAL (kind=8), SAVE :: timpr = 0d0

 !....  PERFORM SELECT TASK

 SELECT CASE (TRIM(itask))

 CASE ('NEW','NSTRA0','NSTRA1','NSTRA2')   !Read Input Data
     !from npo_db INTENT(IN) :: npoin,coord,label,oldlb
     !            INTENT(IN OUT) :: bottom,top,emass
     nn = 0   ! extended check statistics
     ! .. Initializes data base and input element DATA
     CALL cinpu3(maxve,iwrit,npoin,coord,top,bottom)

 CASE ('LUMASS')
   ! compute surface mass
   CALL surms3(emass,coord)

 CASE ('FORCES')  !Performs contact search and computes contact forces
   !from npo_db INTENT(IN) :: npoin,coora,coorb,coort,label,emass
   !            INTENT(IN OUT) :: fcont

   ! ....  compute maximum displacement increment possible
   td = tdtime
   IF( td == 0 ) td = 1d0              !to avoid an error due to initial plastic work
   disma = ABS( MAXVAL(velnp) - MINVAL(velnp) )
   disma = 11.d0*disma * dtcal * SQRT(3d0)  !maximum increment
   IF( ctime < 0d0 )ctime = dtcal      !first step only
   dtime = ctime                       !contact dtime from database
   IF(dtime == 0d0)dtime = dtcal       !computation time

   ! ....  Perform Contact Searching & Compute contact forces
   IF( static .AND. updiv)THEN  !initializes values
     pair => headp
     DO ipair=1,npair                      !for each pair
       IF( pair%press ) pair%presn = 0d0
       pair => pair%next
     END DO
     surtf = 0d0                             !initializes for next period
     timpr = 0d0
   END IF
   CALL celmn3(ttime,dtime,td,disma,coora,emass,fcont,coorb,coort)
   timpr = timpr + dtime               !increase elapsed time
 CASE ('DUMPIN')   !write data to a re-start file
   CALL cdump3(npoin)
 !        WRITE(55,"(i5,3i20)",ERR=9999)(i,nn(i,1:3),i=1,41)  !statistic to debug file

 CASE ('RESTAR')   !read data from a re-start file
   !from npo_db INTENT(IN) :: bottom,top
   CALL crest3(npoin)
   ALLOCATE ( surtf(3,npair) )              !get memory for total forces
   surtf = 0d0                              !initializes

 CASE ('UPDLON')   !Modifies internal numeration
   CALL cupdl3( )

 CASE ('OUTDY1')   ! for History
   !Writes contact forces between pairs
   !Initializes average nodal forces also
   dtime = ctime                       !contact dtime from database
   IF(dtime == 0d0)dtime = dtcal       !computation time
   WRITE(41,ERR=9999) ttime            !control variable
   pair => headp
   DO ipair=1,npair                      !for each pair
     IF(timpr > 0d0 )THEN
       !average contact forces on surfaces
       WRITE(41,ERR=9999) surtf(1:3,ipair)/timpr
       !initializes average contact forces on Nodes
       IF( pair%press ) pair%presn = pair%presn/(timpr+dtime)*dtime
     ELSE
       WRITE(41,ERR=9999) surtf(1:3,ipair)           !average contact forces
     END IF
     pair => pair%next
   END DO
   surtf = 0d0                             !initializes for next period
   timpr = 0d0                             !initializes elapsed time

 CASE ('OUTDY2')   ! for STP
   !Writes contact forces and gaps for nodes
   dtime = ctime                       !contact dtime from database
   IF(dtime == 0d0)dtime = dtcal       !computation time
   auxil = timpr+dtime
   IF( auxil == 0d0 ) auxil = 1d0
   pair => headp
   DO ipair=1,npair                      !for each pair
     IF( pair%press ) WRITE(44,ERR=9999) (pair%presn(i)/auxil,i=1,pair%ncnod)
     IF( pair%wrink ) WRITE(44,ERR=9999) (pair%rssdb(3,i),i=1,pair%ncnod)
     IF( pair%wrink ) WRITE(44,ERR=9999) (pair%mingp(i),i=1,pair%ncnod)
     pair => pair%next
   END DO
   IF( wear ) WRITE(45,ERR=9999) (wwear(i),i=1,npoin)

 CASE ('INITIA')   !Initial penetrations
   ! ....  Perform Contact Searching & Output a report
   CALL celmn3p(ttime,coora,coorb,coort)

 CASE ('WRTSUR')   !Writes contact surfaces for Post-processing
   ! OPEN files for postprocess, and detect surface to use
   CALL prsur3(iwrit)
 END SELECT

 RETURN
 9999 CALL runen2('')

 END SUBROUTINE conta3
