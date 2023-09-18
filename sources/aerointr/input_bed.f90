
SUBROUTINE input_bed
  
  ! Mauro S. Maza - 15/09/2015
  
  ! INPUT Blade Element Distribution
  ! Reads data determining wich elements belong to upper and lower surfaces fo each blade
  
  USE inter_db,     ONLY:   ini_bed, add_bed, & ! Subroutines
                            bed,              & ! Types
                            headbed, tailbed    ! Variables
  
  IMPLICIT NONE
  
  INTEGER             :: j,i,k
  
  TYPE(bed),  POINTER :: newbed
  
  
  OPEN( UNIT = 221, FILE = 'bed.dat', STATUS = 'old', FORM = 'formatted', ACCESS = 'sequential', IOSTAT = j )
    
  IF( j .GT. 0 )THEN
    PRINT *, 'File BED.DAT not found - PROGRAM  s t o p p e d'
  ELSE
    READ( 221 , '(/)' ) ! header
    CALL ini_bed(headbed, tailbed)
    DO i=1,3
      
      ALLOCATE(newbed) ! allocate memory
      CALL add_bed(newbed, headbed, tailbed)
      
      READ( 221 , '(/)' ) ! set name
      READ( 221 , * ) newbed%m_est
      
      READ( 221 , * ) ! upper surface
      READ( 221 , * ) newbed%numuse
      ALLOCATE( newbed%buse(newbed%numuse) )
      DO k=1,newbed%numuse
        READ( 221 , * ) newbed%buse(k)
      ENDDO
      
      READ( 221 , * ) ! lower surface
      READ( 221 , * ) newbed%numlse
      ALLOCATE( newbed%blse(newbed%numlse) )
      DO k=1,newbed%numlse
        READ( 221 , * ) newbed%blse(k)
      ENDDO
      
    ENDDO
    
    CLOSE( UNIT = 221 )
  ENDIF
  
ENDSUBROUTINE input_bed
