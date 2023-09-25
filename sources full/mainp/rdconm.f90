      SUBROUTINE rdconm(actio,ndofn,iwrit)

      !***  READ concentrated masses

      USE c_input
      USE cms_db
      IMPLICIT NONE
      CHARACTER(len=*),INTENT(IN):: actio
      INTEGER (kind=4) :: iwrit,ndofn

      ! Local
      TYPE (cms_nod), POINTER :: cms,head,tail
      INTEGER :: n,rdof

      rdof = MIN(6,ndofn)
      CALL listen('INPDAT')      !read a card

      IF (.NOT.exists('CONCEN')) THEN        !key-word CONCENTRATED_MASSES
        backs = .TRUE.                       !one line back
        RETURN
      END IF
      !Initialize empty list
      CALL ini_cms(head,tail)
      IF (TRIM(actio) == 'NEW' .OR. .NOT.exists('ADD   ')) THEN
        nconm=0
      ELSE IF (exists('ADD   ')) THEN
        DO n=1,nconm
          ALLOCATE (cms)
          cms%node = nodcms(n)
          cms%mass(1:rdof) = cmass(1:rdof,n)
          CALL add_cms (cms, head, tail)
        END DO
      END IF

      IF (ASSOCIATED (nodcms) ) NULLIFY (nodcms)
      IF (ASSOCIATED (cmass) ) NULLIFY (cmass)

      !Loop to read data and add them to the list
      DO
        CALL listen('RDCONM')
        IF (exists('ENDCON')) EXIT

        ALLOCATE (cms)
        cms%node =  INT(param(1))
        cms%mass(1:rdof)= param(2:rdof+1)
        nconm=nconm+1
        CALL add_cms( cms, head, tail )
      END DO

      !Store data in an array
      ALLOCATE (nodcms(nconm), cmass(rdof,nconm))
      CALL store_cms (head,rdof) !transfer to array and release memory

      IF (iwrit == 1) THEN
!        IF(ndime == 2) WRITE(lures,"(/,' Concentrated Masses',
!     &               /,'   Node', 6X,'X',11X,'Y',11X,'RZ')",ERR=9999)
!        IF(ndime == 3)
      WRITE(lures,"(/,' Concentrated Masses', &
     &                               /,'   Node', &
     &              6X,'X',11X,'Y',11X,'z',10x,'RX',10X,'RY',10X,'RZ')",ERR=9999)
        DO n = 1,nconm
          WRITE(lures, '(i8,6e12.4)',ERR=9999) nodcms(n),cmass(1:rdof,n)
        END DO
      END IF

      CALL ini_cms(head,tail)  !nullify pointers

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE rdconm
