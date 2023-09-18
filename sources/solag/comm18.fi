 SUBROUTINE comm18(indic,nelem,nreqs,narch,elsnam,elset,ngaus)

 ! parameters of present set

 IMPLICIT NONE
 CHARACTER(len=*) :: elsnam
 INTEGER (kind=4) :: indic,nelem,nreqs,narch,ngaus
 TYPE (ele18_set) :: elset

 IF(indic == 0) THEN
    elset%sname = elsnam         !set name
    elset%nelem = nelem          !number of elements
    elset%nreqs = nreqs          !number of required elements for output
    elset%narch = narch          !file for output
    elset%ngaus = ngaus          !number of integration points
  ELSE
    elsnam = elset%sname
    nelem  = elset%nelem
    nreqs  = elset%nreqs
    narch  = elset%narch
    ngaus  = elset%ngaus
  END IF

 RETURN
 END SUBROUTINE comm18
