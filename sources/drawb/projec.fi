SUBROUTINE projec (xs, ys, zs, xx1, xx2, xx3, n,  sn, tn, ier)

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER, INTENT (OUT) :: ier
  REAL (kind=8), INTENT (IN) :: xs
  REAL (kind=8), INTENT (IN) :: ys
  REAL (kind=8), INTENT (IN) :: zs
  REAL (kind=8), INTENT (INOUT) :: sn
  REAL (kind=8), INTENT (INOUT) :: tn
  REAL (kind=8), INTENT (IN) :: xx1 (4)
  REAL (kind=8), INTENT (IN) :: xx2 (4)
  REAL (kind=8), INTENT (IN) :: xx3 (4)
  REAL (kind=8), INTENT (INOUT) :: n (3)
  !-----------------------------------------------
  !   l o c a l   v a r i a b l e s

  INTEGER :: i
  REAL (kind=8) :: x13, x42, x31, x24, a1, b1, c1, d1, a2, b2, c2,  &
 &                 d2, a3, b3, c3, d3, xg, yg, zg, e1, e2, e3, a, b,&
 &                 c, r1, r2, odet, pxs, pys, pzs, pxt, pyt, pzt,   &
 &                 sold, told, error, xxx, yyy, zzz, omag
  !-----------------------------------------------------------------

  x13 = xx1 (1) + xx1 (3)
  x42 = xx1 (4) + xx1 (2)
  x31 = xx1 (3) - xx1 (1)
  x24 = xx1 (2) - xx1 (4)
  a1 = .25 * (x13+x42)
  b1 = .25 * (x31+x24)
  c1 = .25 * (x31-x24)
  d1 = .25 * (x13-x42)
  x13 = xx2 (1) + xx2 (3)
  x42 = xx2 (4) + xx2 (2)
  x31 = xx2 (3) - xx2 (1)
  x24 = xx2 (2) - xx2 (4)
  a2 = .25 * (x13+x42)
  b2 = .25 * (x31+x24)
  c2 = .25 * (x31-x24)
  d2 = .25 * (x13-x42)
  x13 = xx3 (1) + xx3 (3)
  x42 = xx3 (4) + xx3 (2)
  x31 = xx3 (3) - xx3 (1)
  x24 = xx3 (2) - xx3 (4)
  a3 = .25 * (x13+x42)
  b3 = .25 * (x31+x24)
  c3 = .25 * (x31-x24)
  d3 = .25 * (x13-x42)
  xg = a1
  yg = a2
  zg = a3
  sn = 0.
  tn = 0.
  e1 = xs - xg
  e2 = ys - yg
  e3 = zs - zg
  a = b1 ** 2 + b2 ** 2 + b3 ** 2
  b = b1 * c1 + b2 * c2 + b3 * c3
  c = c1 ** 2 + c2 ** 2 + c3 ** 2
  r1 = b1 * e1 + b2 * e2 + b3 * e3
  r2 = c1 * e1 + c2 * e2 + c3 * e3
  odet = 1. / (a*c-b**2)
  sn = (r1*c-r2*b) * odet
  tn = (r2*a-r1*b) * odet
  xg = a1 + (b1+d1*tn) * sn + c1 * tn
  yg = a2 + (b2+d2*tn) * sn + c2 * tn
  zg = a3 + (b3+d3*tn) * sn + c3 * tn
  e1 = xs - xg
  e2 = ys - yg
  e3 = zs - zg
  pxs = b1 + d1 * tn
  pys = b2 + d2 * tn
  pzs = b3 + d3 * tn
  pxt = c1 + d1 * sn
  pyt = c2 + d2 * sn
  pzt = c3 + d3 * sn
  a = pxs ** 2 + pys ** 2 + pzs ** 2
  b = pxs * pxt + pys * pyt + pzs * pzt
  c = pxt ** 2 + pyt ** 2 + pzt ** 2
  r1 = pxs * e1 + pys * e2 + pzs * e3
  r2 = pxt * e1 + pyt * e2 + pzt * e3
  odet = 1. / (a*c-b**2)
  sn = sn + (r1*c-r2*b) * odet
  tn = tn + (r2*a-r1*b) * odet
  xg = a1 + (b1+d1*tn) * sn + c1 * tn
  yg = a2 + (b2+d2*tn) * sn + c2 * tn
  zg = a3 + (b3+d3*tn) * sn + c3 * tn
  e1 = xs - xg
  e2 = ys - yg
  e3 = zs - zg
  pxs = b1 + d1 * tn
  pys = b2 + d2 * tn
  pzs = b3 + d3 * tn
  pxt = c1 + d1 * sn
  pyt = c2 + d2 * sn
  pzt = c3 + d3 * sn
  a = pxs ** 2 + pys ** 2 + pzs ** 2
  b = pxs * pxt + pys * pyt + pzs * pzt
  c = pxt ** 2 + pyt ** 2 + pzt ** 2
  r1 = pxs * e1 + pys * e2 + pzs * e3
  r2 = pxt * e1 + pyt * e2 + pzt * e3
  odet = 1. / (a*c-b**2)
  sn = sn + (r1*c-r2*b) * odet
  tn = tn + (r2*a-r1*b) * odet
  xg = a1 + (b1+d1*tn) * sn + c1 * tn
  yg = a2 + (b2+d2*tn) * sn + c2 * tn
  zg = a3 + (b3+d3*tn) * sn + c3 * tn
  ier = 1
  sold = sn
  told = tn
  DO i = 1, 10
     pxs = b1 + d1 * tn
     pys = b2 + d2 * tn
     pzs = b3 + d3 * tn
     pxt = c1 + d1 * sn
     pyt = c2 + d2 * sn
     pzt = c3 + d3 * sn
     e1 = xs - xg
     e2 = ys - yg
     e3 = zs - zg
     a = pxs ** 2 + pys ** 2 + pzs ** 2
     b = pxs *pxt + pys *pyt + pzs * pzt - e1 * d1 - e2 * d2 - e3*d3
     c = pxt ** 2 + pyt ** 2 + pzt ** 2
     r1 = pxs * e1 + pys * e2 + pzs * e3
     r2 = pxt * e1 + pyt * e2 + pzt * e3
     odet = 1. / (a*c-b**2)
     sn = sn + (r1*c-r2*b) * odet
     tn = tn + (r2*a-r1*b) * odet
     xg = a1 + (b1+d1*tn) * sn + c1 * tn
     yg = a2 + (b2+d2*tn) * sn + c2 * tn
     zg = a3 + (b3+d3*tn) * sn + c3 * tn
     error = (sn-sold) ** 2 + (tn-told) ** 2
     IF (error/(1.e-30+sn**2+tn**2) < 2.e-08) GO TO 20
     sold = sn
     told = tn
  END DO
  IF (Abs(sn) > 1.05 .OR. Abs(tn) > 1.05) ier = 2
  GO TO 21
20 CONTINUE
  ier = 0
21 CONTINUE
  IF (Abs(sn) > 1.05 .OR. Abs(tn) > 1.05) ier = 2
  xxx = xg
  yyy = yg
  zzz = zg
  pxs = b1 + d1 * tn
  pys = b2 + d2 * tn
  pzs = b3 + d3 * tn
  pxt = c1 + d1 * sn
  pyt = c2 + d2 * sn
  pzt = c3 + d3 * sn
  n (1) = pys * pzt - pzs * pyt
  n (2) = pzs * pxt - pxs * pzt
  n (3) = pxs * pyt - pys * pxt
  omag = 1. / Sqrt (n(1)**2+n(2)**2+n(3)**2)
  n (1) = omag * n (1)
  n (2) = omag * n (2)
  n (3) = omag * n (3)
  RETURN
END SUBROUTINE projec
