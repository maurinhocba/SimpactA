 SUBROUTINE cspin2(p)

 !.... read contact pair data
 IMPLICIT NONE
 !     arguments
 TYPE (pair2_db), POINTER :: p  !INTENT(OUT)

 !     local variables
 INTEGER (kind=4) j
 CHARACTER(len=mnam):: surpos
 LOGICAL :: found

 !.... read master and slave contact surfaces for each pair
 p%pname = get_name('IPAIR ',found,'!  PAIR NAME:',stype='PAIR')
 CALL add_name(p%pname,4)
 IF (exists('SLAVE ',j)) THEN            !check Key-word
   p%slave = get_name(posin=j,stype='SURF')
   IF( param(j) == 0)THEN                !if no associated number
     surpos  = words(j+2)(1:midn)        !position (top center bottom)
   ELSE                                  !label is a two digit number
     surpos  = words(j+1)(1:midn)        !pos. (top center reversed bottom)
   END IF
   SELECT CASE (TRIM(surpos))
   CASE ('BOTTOM')
     p%slsur = -2
   CASE ('REVERS')
     p%slsur = -1
   CASE ('TOP')
     p%slsur =  2
   CASE DEFAULT
     p%slsur =  1
     surpos = 'CENTER'
   END SELECT
   WRITE(lures,"(9x,'Slave surface :',A,'  Surface posit.:',A)",ERR=9999) &
                 TRIM(p%slave),TRIM(surpos)
 ELSE
   CALL runend ('CSPINP: Slave surface not defined !')
 END IF

 IF (exists('MASTER',j)) THEN           !check Key-word
   p%master = get_name(posin=j,stype='SURF')
   IF( param(j) == 0)THEN               !if no associated number
     surpos  = words(j+2)(1:midn)       !position (top center bottom)
   ELSE                                 !label is a two digit number
     surpos  = words(j+1)(1:midn)       !position (top center bottom)
   END IF
   SELECT CASE (TRIM(surpos))
   CASE ('BOTTOM')
     p%mtsur = -2
   CASE ('REVERS')
     p%mtsur = -1
   CASE ('TOP')
     p%mtsur =  2
   CASE DEFAULT
     p%mtsur =  1
     surpos = 'CENTER'
   END SELECT
   WRITE(lures,"('         Master surface:',a, &
                &'         Surface posit.:',A)",ERR=9999) TRIM(p%master), TRIM(surpos)
 ELSE
   CALL runend ('CSPINP: Master surface not defined!')
 END IF

 p%auto = (p%master == p%slave) .AND. (p%mtsur == p%slsur)  !self contact

 !.... read contact pair properties data

 p%indcon =getint('INDCON',0,'!Contact forces flag ..............')
 IF( p%auto .AND. p%indcon /= 0 )THEN
   p%indcon = 0
   WRITE(lures,"(15x,'WARNING For SELF Contact Forces must be applied',/, &
     &           15x,' on BOTH surfaces, CHANGED')",ERR=9999)
 END IF
 SELECT CASE (p%indcon)
 CASE (0)
   WRITE(lures,"(15x,'Contact Forces applied on BOTH surfaces')",ERR=9999)
 CASE (1)
   WRITE(lures,"(15x,'Contact Forces applied on SLAVE surface ONLY')",ERR=9999)
 CASE (2)
   WRITE(lures,"(15x,'Contact Forces applied on MASTER surface ONLY')",ERR=9999)
 CASE DEFAULT
   p%indcon = 0
   WRITE(lures,"(15x,'Contact Forces applied on BOTH surfaces')",ERR=9999)
 END SELECT

 p%npenal =getrea('NPENAL',0d0,'!Normal penalty coeff .............')
 p%tpenal =getrea('TPENAL',0d0,'!Tangential penalty coeff..........')
 p%static =getrea('STATIC',0d0,'!Static Friction coefficient ......')
 p%kinet  =getrea('KINET ',p%static,            &
                  ' Kinetic Friction coefficient .....')
 p%cutoff=-getrea('CUTOFF',1d9,' CUT OFF value for penetration.....')
 IF( p%cutoff == 0d0 ) p%cutoff = -1d9
 p%gapinc=-getrea('GAPINC',-p%cutoff/20,        &
                  ' Maximum Incremental Gap per Step..')
 IF( p%gapinc >= 0d0 ) p%gapinc = p%cutoff/20
 IF( therm )THEN
   p%hcont  =getrea('HCONT ',0d0,' H factor for heat conduction......')
   p%econt = getrea('EHCON ',1d0,' Exponent for heat conduction......')
   p%hv    = getrea('VICKH ',1d0,' Vickers  hardness ................')
   p%rec   = getrea('EFFUS ',5d-1,' Relative Effusivity (Slave).......')
 ELSE
   p%hcont = 0d0
   p%econt = 0d0
   p%hv    = 0d0
   p%rec   = 0d0
 END IF
 p%start  =getrea('START ',0d0,' Activation time ..................')
 p%end    =getrea('END   ',1d9,' Deactivation time  ...............')
 IF( p%auto )THEN
   p%freq   =getint('FREQ  ',1000,' Frequency of global search .......')
 ELSE
   p%freq   =getint('FREQ  ',10000,' Frequency of global search .......')
 END IF
 p%cpress = exists('CPRESS')
 p%press = exists('PRESS ') .OR. p%cpress
 IF( p%press ) WRITE(lures,"(15x,'Contact pressure will be computed')",ERR=9999)

 p%wrink = exists('WRINKL')
 IF( p%wrink ) WRITE(lures,"(15x,'Wrinkle control will be performed')",ERR=9999)

 p%prev = .FALSE.  !not existence in previous step
 WRITE(lures,"(//)",ERR=9999)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE cspin2
