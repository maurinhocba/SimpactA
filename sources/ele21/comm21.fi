SUBROUTINE comm21(indic,nelem,nreqs,narch,nnode,ngaus,  &
                  elsnam,elset)

  IMPLICIT NONE

  CHARACTER(len=*) :: elsnam
  INTEGER (kind=4) :: indic,nelem,nreqs,nnode,ngaus,narch
  TYPE (ele21_set) :: elset

  IF(indic == 0) THEN
     elset%sname  = elsnam
     elset%nelem  = nelem
     elset%nreqs  = nreqs
     elset%ngaus  = ngaus
     elset%nnode  = nnode
     elset%narch  = narch
  ELSE
     elsnam = elset%sname
     nelem  = elset%nelem
     nreqs  = elset%nreqs
     ngaus  = elset%ngaus
     nnode  = elset%nnode
     narch  = elset%narch
  END IF

END SUBROUTINE comm21
