 SUBROUTINE cspin1(p)

 !.... read contact pair data
 IMPLICIT NONE
 !     arguments
 TYPE (pair1_db), POINTER :: p  !INTENT(OUT)

 !     local variables
 INTEGER (kind=4) j
 LOGICAL :: found

 !.... read master and slave contact surfaces for each pair
 p%pname = get_name('IPAIR ',found,'!  PAIR NAME:',stype='PAIR')
 CALL add_name(p%pname,4)
 IF (exists('SLAVE ',j)) THEN            !check Key-word
   p%slave = get_name(posin=j,stype='SURF')
   WRITE(lures,"(9x,'Slave surface :',A)",ERR=9999)  TRIM(p%slave)
 ELSE
   CALL runend ('CSPINP: Slave surface not defined !')
 END IF

 IF (exists('MASTER',j)) THEN           !check Key-word
   p%master = get_name(posin=j,stype='SURF')
   WRITE(lures,"('         Master surface:',a)",ERR=9999) TRIM(p%master)
 ELSE
   CALL runend ('CSPINP: Master surface not defined!')
 END IF

 !.... read contact pair properties data

 p%indcon =getint('INDCON',0,'!Contact forces flag ..............')
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
 p%kinet  =getrea('KINET ',p%static,' Kinetic Friction coefficient .....')
 p%start  =getrea('START ',0d0,' Activation time ..................')
 p%end    =getrea('END   ',1d9,' Deactivation time  ...............')
 p%press = exists('PRESS ')
 IF( p%press ) WRITE(lures,"(15x,'Contact pressure will be computed')",ERR=9999)
 p%prev = .FALSE.  !not existence in previous step
 WRITE(lures,"(//)",ERR=9999)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE cspin1
