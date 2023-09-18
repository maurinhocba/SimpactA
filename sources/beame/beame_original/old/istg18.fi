 SUBROUTINE istg18(ng,nstre,stres,dstra,epdef,sedef,props,istop,ielem)
 !**********************************************************************
 !
 !      incremental stresses at gauss points for beams (model 1)
 !
 !**********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) ng,nstre,istop,ielem
 REAL    (kind=8) stres(nstre,ng),dstra(nstre,ng),epdef(ng), &
                  sedef(ng),props(*),ws(36)
 INTEGER (kind=4) g,i
 REAL    (kind=8) def

 DO g=1,ng
   sedef(g) = sedef(g) + dstra(1,g)
   def   = 1d0/sedef(g)
   ws(1) = props(1)*def
   ws(2) = props(2)*def
   ws(3) = props(3)*def
   ws(4) = props(4)*def**2
   ws(5) = props(5)*def**2
   ws(6) = props(6)*def**2
   DO i=1,nstre
     stres(i,g) = stres(i,g) + ws(i)*dstra(i,g)
   END DO
   !       modify the plasticity constants to INCLUDE changes in section)
   IF(props(7) > 0) THEN
     ws(7) = props(7)*def
     ws(8) = props(8)*def
     ws(9) = props(9)*def
     def   = def**(1.5d0)
     ws(10) = props(10)*def
     ws(11) = props(11)*def
     ws(12) = props(12)*def
     def = props(17)**2
     CALL rare18(stres(1,g),ws,epdef(g),props(13),props(14), &
                 props(16),props(17),ws(15),ws(31),istop,ielem)
     IF(istop == 1) RETURN
   END IF
 END DO

 END SUBROUTINE istg18
