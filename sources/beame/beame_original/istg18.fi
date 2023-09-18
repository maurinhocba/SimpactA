 SUBROUTINE istg18(ng,nstre,stres,dstra,epdef,sedef,props,istop,ielem,plast)
 !**********************************************************************
 !
 !      incremental stresses at gauss points for beams (model 1)
 !
 !**********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) ng,nstre,istop,ielem
 REAL    (kind=8) stres(nstre,ng),dstra(nstre,ng),epdef(ng), &
                  sedef(ng),props(:),ws(36)
 LOGICAL :: plast
 INTEGER (kind=4) g,i
 REAL    (kind=8) def

 DO g=1,ng
   IF( plast )THEN
     sedef(g) = sedef(g) + dstra(1,g)
   ELSE
     sedef(g) = 1d0 + dstra(1,g)
   END IF
   def   = 1d0/sedef(g)
   ws(1) = props(5)*def
   ws(2) = props(6)*def
   ws(3) = props(7)*def
   ws(4) = props(8)*def**2
   ws(5) = props(9)*def**2
   ws(6) = props(10)*def**2
   !       modify the plasticity constants to INCLUDE changes in section)
   IF( plast )THEN
     DO i=1,nstre
       stres(i,g) = stres(i,g) + ws(i)*dstra(i,g)
     END DO
     ws(7) = props(11)*def
     ws(8) = props(12)*def
     ws(9) = props(13)*def
     def   = def**(1.5d0)
     ws(10) = props(14)*def
     ws(11) = props(15)*def
     ws(12) = props(16)*def
     def = props(21)**2
     CALL rare18(stres(1,g),ws,epdef(g),props(17),props(18), &
                 props(20),props(21),ws(15),ws(31),istop,ielem)
     IF(istop == 1) RETURN
   ELSE
     DO i=1,nstre
       stres(i,g) = ws(i)*dstra(i,g)
     END DO
   END IF
 END DO

 END SUBROUTINE istg18
