 SUBROUTINE comm11 (indic, nelem, nstre, nbn, ngaus, &
                    nreqs, narch, elsnam, elset)

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 INTEGER (kind=4) :: indic, nelem, nstre, nbn, &
&                    nreqs, narch, ngaus
 TYPE (ele11_set) :: elset


 IF(indic == 0) THEN
    elset%sname  = elsnam
    elset%nelem  = nelem
    elset%nstre  = nstre
    elset%nbn    = nbn
    elset%nreqs  = nreqs
    elset%narch  = narch
    elset%ngaus  = ngaus
 ELSE
    elsnam = elset%sname
    nelem = elset%nelem
    nstre = elset%nstre
    nbn   = elset%nbn
    nreqs = elset%nreqs
    narch = elset%narch
    ngaus = elset%ngaus
 END IF

 END SUBROUTINE comm11
