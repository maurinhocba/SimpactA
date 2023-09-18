 SUBROUTINE commv6 (indic, nelem, nreqs, narch, elsnam, elset)

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 INTEGER (kind=4) :: indic,nelem, nreqs, narch
 TYPE (ele06_set) :: elset

 IF(indic == 0) THEN
    elset%sname  = elsnam
    elset%nelem  = nelem
    elset%nreqs  = nreqs
    elset%narch  = narch
 ELSE
    nelem  = elset%nelem
    elsnam = elset%sname
    nreqs   = elset%nreqs
    narch   = elset%narch
 END IF

 END SUBROUTINE commv6
