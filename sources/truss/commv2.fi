 SUBROUTINE commv2(indic,nelem,nreqs,narch,elsnam,elset)

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 INTEGER (kind=4) :: indic,nelem,nreqs,narch
 TYPE (ele02_set) :: elset

 IF(indic == 0) THEN
    elset%sname = elsnam
    elset%nelem = nelem
    elset%nreqs = nreqs
    elset%narch = narch
  ELSE
    nelem = elset%nelem
    nreqs = elset%nreqs
    narch = elset%narch
    elsnam = elset%sname
  END IF

 END SUBROUTINE commv2
