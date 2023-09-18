 SUBROUTINE sc_spot ( sec)
 !            -------------------------------------------------
 !            input properties for a point stiffness and damper
 !
 IMPLICIT NONE
 TYPE (section), POINTER :: sec
 !local variales
 LOGICAL :: found
 INTEGER (kind=4) matno,nc,i,k
 INTEGER (kind=4) :: matsa = 20
 CHARACTER (len=mnam) :: cname
 TYPE (curve), POINTER :: c

 CALL listen('SC_SPO')

 sec%secty = 2  ! spot
 matno    = getint('NUMAT ',   1,'!Associated Material ..............')
 ! search for associated material
 CALL mat_search (matno,found,sec%mtbas)
 IF(.NOT.found) CALL runend ('SC_SPO: NOT FOUND ASSOCIATED MATER ')

 IF( sec%mtbas%mtype /= matsa)  &
      CALL runend('SC_SPO: Invalid material association')

 sec%mabas = matno  !material label
 sec%secdef = (/ 6,8,0,0,0 /)
 ALLOCATE (sec%iprop(sec%secdef(1)), sec%rprop(sec%secdef(2)) )
 sec%iprop = 0
 sec%rprop = 0d0

 nc = sec%mtbas%matdef(12)      ! number of curves

 ! read curve for stiffness
 sec%rprop(1) = getrea('ASTIFF',0d0,' Axial Stiffness ..................')
 IF( sec%rprop(1) <= 0d0 .AND. exists('ASTIFF') )THEN
   cname = get_name('ASTIFF') !read curve name
   CALL find_cr(sec%mtbas%chead,c,cname,nc,k)
   sec%iprop(1) = k
   i = 1
   sec%rprop(1) = inte_cr(c%val,c%np,0d0,i)
 END IF
 sec%rprop(2) = getrea('RSTIFF',0d0,' Rotational Stiffness..............')
 IF( sec%rprop(2) <= 0d0 .AND. exists('RSTIFF') )THEN
   cname = get_name('RSTIFF') !read curve name
   CALL find_cr(sec%mtbas%chead,c,cname,nc,k)
   sec%iprop(2) = k
   i = 1
   sec%rprop(2) = inte_cr(c%val,c%np,0d0,i)
 END IF
 IF( exists('XYZ') )THEN
   sec%secdef(3) = 0           !default
 ELSE IF( exists('XY') )THEN
   sec%secdef(3) = 4
 ELSE IF( exists('XZ') )THEN
   sec%secdef(3) = 5
 ELSE IF( exists('YZ') )THEN
   sec%secdef(3) = 6
 ELSE IF( exists('X') )THEN
   sec%secdef(3) = 1           !
 ELSE IF( exists('Y') )THEN
   sec%secdef(3) = 2
 ELSE IF( exists('Z') )THEN
   sec%secdef(3) = 3
 END IF
 ! read curve for damping
 sec%rprop(3) = getrea('ADAMPI',0d0,' Axial Damping ....................')
 IF( sec%rprop(3) <= 0d0 .AND. exists('ADAMPI') )THEN
   cname = get_name('ADAMPI') !read curve name
   CALL find_cr(sec%mtbas%chead,c,cname,nc,k)
   sec%iprop(3) = k
   i = 1
   sec%rprop(3) = inte_cr(c%val,c%np,0d0,i)
 END IF
 sec%rprop(4) = getrea('RDAMPI',0d0,' Rotational Damping ...............')
 IF( sec%rprop(4) <= 0d0 .AND. exists('RDAMPI') )THEN
   cname = get_name('RDAMPI') !read curve name
   CALL find_cr(sec%mtbas%chead,c,cname,nc,k)
   sec%iprop(4) = k
   i = 1
   sec%rprop(4) = inte_cr(c%val,c%np,0d0,i)
 END IF

 ! read mass for translation and rotation
 sec%rprop(5) = getrea('TMASS ',0d0,' Translational Mass................')
 IF( sec%rprop(5) <= 0d0 .AND. exists('TMASS ') )THEN
   cname = get_name('TMASS ') !read curve name
   CALL find_cr(sec%mtbas%chead,c,cname,nc,k)
   sec%iprop(5) = k
   i = 1
   sec%rprop(5) = inte_cr(c%val,c%np,0d0,i)
 END IF
 sec%rprop(6) = getrea('RMASS ',0d0,' Rotational Mass ..................')
 IF( sec%rprop(6) <= 0d0 .AND. exists('RMASS ') )THEN
   cname = get_name('RMASS ') !read curve name
   CALL find_cr(sec%mtbas%chead,c,cname,nc,k)
   sec%iprop(6) = k
   i = 1
   sec%rprop(6) = inte_cr(c%val,c%np,0d0,i)
 END IF

 ! read values for 6-node spots
 sec%rprop(7) = getrea('AREA  ',0d0,' Transversal Area .................')
 sec%rprop(8) = getrea('NPENAL',1d-2,' Normal Penalty ...................')

 RETURN

 END SUBROUTINE sc_spot
