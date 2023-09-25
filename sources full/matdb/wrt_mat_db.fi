 SUBROUTINE wrt_mat_db ( )
 ! write material data necessary for post-process
 ! presently only FLC curves are written

 INTEGER (kind=4), ALLOCATABLE :: flclbl(:),flcs(:)
 TYPE (section), POINTER :: sec         !pointer to a section data
 INTEGER (kind=4) :: i,nfl,j
 ALLOCATE( flclbl(nusect), flcs(nflc) )

 flcs = -1                       !flag
 nfl = 0                         !initializes
 OUTER : DO i=1,nusect
   sec => psecs(i)%p             !point to section
   IF( sec%mabas < 0 )CYCLE OUTER   !non-homogeneuos material
   SELECT CASE (sec%secty)    !section type
   CASE (0,11,12)                !solid or shells
     flclbl(i) = sec%iprop(3)
   CASE DEFAULT                  !rest of sections
     flclbl(i) = 0
   END SELECT
   DO j=1,nfl
     IF( flcs(j) == flclbl(i)) CYCLE OUTER
   END DO
   nfl = nfl + 1
   flcs(nfl) = flclbl(i)
 END DO  OUTER
 WRITE(17,ERR=9999) flclbl,nfl
 DO i=1,nfl
   CALL wrflc(flcs(i),17)      !Write FLC curve on file 17
 END DO

 DEALLOCATE(flclbl,flcs)
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE wrt_mat_db
