 SUBROUTINE conta2(itask,dtcal,ttime,iwrit,velNP,maxve)

  !     main contac routine (ALGORITHM 2)

 USE ctrl_db, ONLY: bottom, top, tdtime
 USE npo_db, ONLY : emass, fcont, coord, coorb, coort
 USE cont2_db
 IMPLICIT NONE

   !        Dummy arguments
                                 !task to perform
   CHARACTER(len=*),INTENT(IN):: itask
   INTEGER(kind=4),INTENT(IN):: iwrit,maxve
   REAL(kind=8),INTENT(IN),OPTIONAL:: dtcal,ttime,velnp(:,:)
   !        Local variables
   TYPE (pair2_db), POINTER :: pair
   INTEGER (kind=4) :: ipair,i
   REAL (kind=8) :: disma, dtime, auxil, td
   REAL (kind=8), SAVE :: timpr = 0d0

   !....  PERFORM SELECT TASK

   SELECT CASE (TRIM(itask))

   CASE ('NEW','NSTRA0','NSTRA1','NSTRA2')   !Read Input Data
       !from npo_db INTENT(IN) :: npoin,coord,label
       !            INTENT(IN OUT) :: bottom,top,emass
       ! .. Initializes data base and input element DATA
       CALL cinpu2(maxve,iwrit,npoin,coord,top,bottom)

   CASE ('LUMASS')
     !  compute mass
     CALL surms2(emass,coord)

   CASE ('FORCES')  !Performs contact search and computes contact forces
     !from npo_db INTENT(IN) :: npoin,coora,coorb,coort,label,emass
     !            INTENT(IN OUT) :: fcont
     ! ....  compute maximum displacement increment possible
     td = tdtime
     IF( td == 0 ) td = 1d0              !to avoid an error due to initial plastic work
     disma = ABS( MAXVAL(velnp) - MINVAL(velnp) )
     disma = 11.d0*disma * dtcal * SQRT(2d0)  !maximum increment
     IF( ctime < 0d0 )ctime = dtcal      !first step only
     dtime = ctime                       !contact dtime from database
     IF(dtime == 0d0)dtime = dtcal       !computation time

     ! ....  Perform Contact Searching & Compute contact forces
     CALL celmn2(ttime,dtime,td,disma,coora,emass,fcont,coorb,coort)
     timpr = timpr + dtime               !increase elapsed time
   CASE ('DUMPIN')   !write data to a re-start file
     CALL cdump2(npoin)

   CASE ('RESTAR')   !read data from a re-start file
     !from npo_db INTENT(IN) :: bottom,top
     CALL crest2(npoin)
     ALLOCATE ( surtf(2,npair) )              !get memory for total forces
     surtf = 0d0                              !initializes

   CASE ('UPDLON')   !Modifies internal numeration
     CALL cupdl2( )

   CASE ('OUTDY1')   !Writes contact forces between pairs
     dtime = ctime                       !contact dtime from database
     IF(dtime == 0d0)dtime = dtcal       !computation time
     WRITE(41,ERR=9999) ttime            !control variable
     pair => headp
     DO ipair=1,npair                      !for each pair
       IF(timpr > 0d0 )THEN
         WRITE(41,ERR=9999) surtf(1:2,ipair)/timpr     !average contact forces
       ELSE
         WRITE(41,ERR=9999) surtf(1:2,ipair)           !average contact forces
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
       IF( pair%wrink ) WRITE(44,ERR=9999) (pair%rssdb(2,i),i=1,pair%ncnod)
       IF( pair%wrink ) WRITE(44,ERR=9999) (pair%mingp(i),i=1,pair%ncnod)
       pair => pair%next
     END DO
     IF( wear ) WRITE(45,ERR=9999) (wwear(i),i=1,npoin)
   CASE ('INITIA')   !Control Initial penetrations
     ! ....  Perform Contact Searching & Output a report
     !print *,' va a entra a celmn2p'
     CALL celmn2p(ttime,coora,coorb,coort)
     !print *,' salio de INITIAL CONTACT'

   CASE ('WRTSUR')   !Writes contact surfaces for Post-processing
     ! OPEN files for postprocess, detect surface to use
     CALL prsur2(iwrit)
   !CASE ('RMMODF')  !Correct body boundary domain during remeshing
   !  CALL rmcont2(sname)
   END SELECT

 RETURN
 9999 CALL runen2('')

 END SUBROUTINE conta2
