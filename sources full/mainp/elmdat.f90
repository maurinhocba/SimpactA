SUBROUTINE elmdat(task,nelem,elsnam,itype)

! read element sets

  USE ctrl_db, ONLY: dtuser, ndime, ndofn, npoin, neulr
  USE outp_db, ONLY: iwrit
  USE esets_db, ONLY: nelms,nel
  USE npo_db
  IMPLICIT NONE

  CHARACTER(len=*),INTENT(IN) :: elsnam
  CHARACTER(len=*),INTENT(IN) :: task
  INTEGER (kind=4),INTENT(IN) :: itype
  INTEGER (kind=4),INTENT(IN OUT) :: nelem

  INTERFACE
    INCLUDE 'inpda1.h'
    INCLUDE 'inpda2.h'
    INCLUDE 'inpd04.h'
    INCLUDE 'inpd05.h'
    INCLUDE 'inpda6.h'
    INCLUDE 'inpda7.h'
    INCLUDE 'inpda8.h'
    INCLUDE 'inpda9.h'
    INCLUDE 'inpd10.h'
    INCLUDE 'inpd11.h'
    INCLUDE 'inpd13.h'
    INCLUDE 'inpd14.h'
    INCLUDE 'inpd15.h'
    INCLUDE 'inpd16.h'
    INCLUDE 'inpd17.h'
    INCLUDE 'inpd18.h'
    INCLUDE 'inpd19.h'
    INCLUDE 'inpd20.h'
    INCLUDE 'inpd21.h'
    INCLUDE 'inpd22.h'
    INCLUDE 'inpd24.h'
    INCLUDE 'inpd25.h'
  END INTERFACE

  SELECT CASE (itype)
  CASE (1)
    CALL inpda1(task,ndime,neulr,nelem,iwrit,elsnam,nelms(1))
  CASE (2)
    CALL inpda2(task,nelem,iwrit,elsnam,nelms(2))
  CASE (4)
    CALL inpd04(task,iwrit,elsnam,nelms(4))
  CASE (5)
    CALL inpd05(task,iwrit,elsnam,nelms(5))
  CASE (6)
    CALL inpda6(task,nelem,eule0,euler,coord,iwrit,elsnam,nelms(6))
  CASE (7)
    CALL inpda7(task,nelem,eule0,euler,coord,iwrit,elsnam,nelms(7))
  CASE (8)
    CALL inpda8(task,nel(8),eule0,euler,coord,iwrit,elsnam,nelms(8))
  CASE (9)
    CALL inpda9(task,nelem,eule0,euler,coord,iwrit,elsnam,nelms(9))
  CASE (10)
    CALL inpd10(task,ndime,nelem,iwrit,elsnam,nelms(10))
  CASE (11)
    CALL inpd11(task,iwrit,elsnam,nelms(11))
  CASE (13)
    CALL inpd13(task,iwrit,elsnam,nelms(13))
  CASE (14)
    CALL inpd14(task,iwrit,elsnam,nelms(14))
  CASE (15)
    CALL inpd15(task,iwrit,elsnam,nelms(15))
  CASE (16)
    CALL inpd16(task,iwrit,elsnam,nelms(16))
  CASE (17)
    CALL inpd17(task,nel(17),iwrit,elsnam,nelms(17))
  CASE (18)
    CALL inpd18(task,nel(18),iwrit,elsnam,nelms(18))
  CASE (19)
    CALL inpd19(task,iwrit,elsnam,nelms(19))
  CASE (20)
    CALL inpd20(task,iwrit,elsnam,nelms(20))
  CASE (21)
    CALL inpd21 (task, nelem, nelms(21), iwrit, elsnam)
  CASE (22)
    CALL inpd22 (task, nelem, nelms(22), iwrit, elsnam)
  CASE (24)
    CALL inpd24(task,iwrit,elsnam,nelms(24))
  CASE (25)
    CALL inpd25(task,iwrit,elsnam,nelms(25))
  CASE DEFAULT
    CALL runend('ELMDAT: ELEMENT_TYPE NOT EXISTENT  ')
  END SELECT

  RETURN
END SUBROUTINE elmdat
