SUBROUTINE trian (nnode, lnod, l, n)

  ! triangularize old mesh

  IMPLICIT NONE

  INTEGER (kind=4) :: n, nnode, lnod(nnode), l(3,2)

  IF (nnode == 4) THEN
    !  4 node --> 2*3 node.
    l(1,1) = lnod(1)
    l(2,1) = lnod(2)
    l(3,1) = lnod(3)
    l(1,2) = lnod(1)
    l(2,2) = lnod(3)
    l(3,2) = lnod(4)
    n = 2
  ELSE IF (nnode == 3) THEN
    !  3 node --> 1*3 node.
    l(1,1) = lnod(1)
    l(2,1) = lnod(2)
    l(3,1) = lnod(3)
    n = 1
  END IF

  RETURN
END
