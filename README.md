#LOCAL REAL Lc
#LOCAL REAL Tc
#LOCAL REAL wt
#LOCAL INTEGER mc
#LOCAL INTEGER nc
#LOCAL REAL deltau

#LOCAL REAL udc_sum1
#LOCAL REAL udc_PI1
#LOCAL REAL udc_sum2
#LOCAL REAL udc_PI2
#LOCAL REAL udc_sum3
#LOCAL REAL udc_PI3

#LOCAL REAL udc_suma(12),udc_PIa(12),udc_sumb(12),udc_PIb(12),da_ac,udc_sumc(12),udc_PIc(12)

#LOCAL REAL db_ac
#LOCAL REAL dc_ac
#LOCAL REAL datmp
#LOCAL REAL dbtmp
#LOCAL REAL dctmp
#LOCAL REAL dmax
#LOCAL REAL dmid
#LOCAL REAL dmin
#LOCAL REAL d0
#LOCAL REAL ims
#LOCAL REAL ima,iqa
#LOCAL REAL imb,iqb
#LOCAL REAL imc,iqc
#LOCAL REAL uabuf
#LOCAL REAL ubbuf
#LOCAL REAL ucbuf
#LOCAL REAL iabuf
#LOCAL REAL ibbuf
#LOCAL REAL icbuf
#LOCAL REAL deltaia
#LOCAL REAL deltaib
#LOCAL REAL deltaic

#LOCAL INTEGER M2,N2,J2
#LOCAL REAL A2(12),B2(12),K2,S2(12),UP2(12),UP3(12),SUM2
#LOCAL REAL A4(12),B4(12),S4(12),SUM4,UP4(12)
#LOCAL REAL A3(12),B3(12),S3(12),SUM3
#LOCAL REAL SUM5,SUM6,SUM7

#LOCAL REAL dua
#LOCAL REAL dub
#LOCAL REAL duc

#LOCAL REAL usab
#LOCAL REAL usbb
#LOCAL REAL uscb

#LOCAL REAL PI

! ???????
DO J2=1,$dim
  S2(J2) = 0.0
  S3(J2) = 0.0
  S4(J2) = 0.0
ENDDO

Lc = $Lm
Tc = 1.0/$fc
mc = floor(TIME/Tc)
M2 = floor(TIME/$DT)
K2 = $T / $DT
DO J2=1,$dim
  A2(J2) = $IN(J2)
  B2(J2) = DELAY($T,64,$IN(J2))

  A3(J2) = $INB(J2)
  B3(J2) = DELAY($T,64,$INB(J2))

  A4(J2) = $INC(J2)
  B4(J2) = DELAY($T,64,$INC(J2))
ENDDO

if (mc > nc) then

  IF (mc .LE. K2) THEN
    DO J2=1,$dim
      S2(J2) = S2(J2) + A2(J2)
      S3(J2) = S3(J2) + A3(J2)
      S4(J2) = S4(J2) + A4(J2)
    ENDDO
  ELSE
    SUM2 = 0.0
    SUM3 = 0.0 
    SUM4 = 0.0 
    DO J2=1,$dim
      S2(J2) = S2(J2) + A2(J2) - B2(J2)
      S3(J2) = S3(J2) + A3(J2) - B3(J2)
      S4(J2) = S4(J2) + A4(J2) - B4(J2)

      UP2(J2) = S2(J2) / K2
      SUM2 = SUM2 + UP2(J2)

      UP3(J2) = S3(J2) / K2
      SUM3 = SUM3 + UP3(J2)

      UP4(J2) = S4(J2) / K2
      SUM4 = SUM4 + UP4(J2)
    ENDDO
  ENDIF

  SUM5 = 0.0
  SUM6 = 0.0
  SUM7 = 0.0
  DO J2=1,$dim
    SUM5 = SUM5 + $IN(J2)
    SUM6 = SUM6 + $INB(J2)
    SUM7 = SUM7 + $INC(J2)
  ENDDO

  $udca = SUM5
  $udcb = SUM6
  $udcc = SUM7

  usab = $ua - 0.10181 * ($ub - $uc)
  usbb = $ub - 0.10181 * ($uc - $ua)
  uscb = $uc - 0.10181 * ($ua - $ub)

  $x = usab

  dua = $udcref - $izpa / 12
  if (dua > 0.15) then
    dua = 0.15
  elseif (dua < -0.15) then
    dua = -0.15
  endif
  udc_sum1 = udc_sum1 + $Ki_udc * dua * Tc
  if (udc_sum1 > $Integ_udc_max) then
    udc_sum1 = $Integ_udc_max
  elseif (udc_sum1 < -$Integ_udc_max) then
    udc_sum1 = -$Integ_udc_max
  endif
  udc_PI1 = $Kp_udc * dua + udc_sum1
  if (udc_PI1 > $PI_udc_max) then
    udc_PI1 = $PI_udc_max
  elseif (udc_PI1 < -$PI_udc_max) then
    udc_PI1 = -$PI_udc_max
  endif
  if (TIME < 0.812) then
    udc_sum1 = 0.0
    udc_PI1 = 0.0
  endif

  ima = udc_PI1 * usab / 2.2

  dub = $udcref - $izpb / 12
  if (dub > 0.15) then
    dub = 0.15
  elseif (dub < -0.15) then
    dub = -0.15
  endif
  udc_sum2 = udc_sum2 + $Ki_udc * dub * Tc
  if (udc_sum2 > $Integ_udc_max) then
    udc_sum2 = $Integ_udc_max
  elseif (udc_sum2 < -$Integ_udc_max) then
    udc_sum2 = -$Integ_udc_max
  endif
  udc_PI2 = $Kp_udc * dub + udc_sum2
  if (udc_PI2 > $PI_udc_max) then
    udc_PI2 = $PI_udc_max
  elseif (udc_PI2 < -$PI_udc_max) then
    udc_PI2 = -$PI_udc_max
  endif
  if (TIME < 0.812) then
    udc_sum2 = 0.0
    udc_PI2 = 0.0
  endif
  imb = udc_PI2 * usbb / 2.2

  duc = $udcref - $izpc / 12
  if (duc > 0.15) then
    duc = 0.15
  elseif (duc < -0.15) then
    duc = -0.15
  endif
  udc_sum3 = udc_sum3 + $Ki_udc * duc * Tc
  if (udc_sum3 > $Integ_udc_max) then
    udc_sum3 = $Integ_udc_max
  elseif (udc_sum3 < -$Integ_udc_max) then
    udc_sum3 = -$Integ_udc_max
  endif
  udc_PI3 = $Kp_udc * duc + udc_sum3
  if (udc_PI3 > $PI_udc_max) then
    udc_PI3 = $PI_udc_max
  elseif (udc_PI3 < -$PI_udc_max) then
    udc_PI3 = -$PI_udc_max
  endif
  if (TIME < 0.812) then
    udc_sum3 = 0.0
    udc_PI3 = 0.0
  endif
  imc = udc_PI3 * uscb / 2.2

  PI = 3.141592653589793
  wt = 2.0 * PI * 50.0 * TIME

  iqa = 0 * (usbb - uscb) / 1000
  iqb = 0 * (uscb - usab) / 1000
  iqc = 0 * (usab - usbb) / 1000

  deltaia = 0.45 * (-ima + $ira + iqa - $ia)
  deltaib = 0.45 * (-imb + $irb + iqb - $ib)
  deltaic = 0.45 * (-imc + $irc + iqc - $ic)

  da_ac = (deltaia * Lc / Tc + usab / 2.2) / SUM5
  db_ac = (deltaib * Lc / Tc + usbb / 2.2) / SUM6
  dc_ac = (deltaic * Lc / Tc + uscb / 2.2) / SUM7
if (da_ac > db_ac) then
    dmax = da_ac
    dmin = db_ac
  else
    dmax = db_ac
    dmin = da_ac
  endif
  if (dc_ac > dmax) then
    dmax = dc_ac
  endif
  if (dc_ac < dmin) then
    dmin = dc_ac
  endif

  d0 = -0.5 * (dmax + dmin)

  datmp = da_ac + d0
  dbtmp = db_ac + d0
  dctmp = dc_ac + d0

endif
#LOCAL REAL z1,z2,z3
#LOCAL REAL y1,y2,y3
#LOCAL REAL T
#LOCAL INTEGER M,N

      T=1.0/$fs
      M=FLOOR(TIME/T)
      IF(M>N) THEN
      z3=z2
      z2=z1
      z1=$eaz
      y3=y2
      y2=y1
      y1=-$a2*y2-$a3*y3+$b1*z1+$b2*z2+$b3*z3
      $eaf=y1
      N=M
      ENDIF
#LOCAL REAL z11,z12,z13
#LOCAL REAL y11,y12,y13
#LOCAL REAL T1
#LOCAL INTEGER M1,N1

      T1=1.0/$fs
      M1=FLOOR(TIME/T1)
      IF(M1>N1) THEN
      z13=z12
      z12=z11
      z11=$ebz
      y13=y12
      y12=y11
      y11=-$a2*y12-$a3*y13+$b1*z11+$b2*z12+$b3*z13
      $ebf=y11
      N1=M1
      ENDIF
#LOCAL REAL z16,z26,z36,z46
#LOCAL REAL y16,y26,y36,y46
#LOCAL REAL T6
#LOCAL INTEGER M6,N6

      T6=1.0/$fs
      M6=FLOOR(TIME/T6)
      IF(M6>N6) THEN
      z46=z36
      z36=z26
      z26=z16
      z16=$eaaz
      y46=y36
      y36=y26
      y26=y16
      y16=-$a2*y26-$a3*y36-$a4*y46+$b1*z16+$b2*z26+$b3*z36+$b4*z46
      $eaaf=y16
      N6=M6
      ENDIF
#LOCAL REAL z17,z27,z37,z47
#LOCAL REAL y17,y27,y37,y47
#LOCAL REAL T7
#LOCAL INTEGER M7,N7

      T7=1.0/$fs
      M7=FLOOR(TIME/T7)
      IF(M7>N7) THEN
      z47=z37
      z37=z27
      z27=z17
      z17=$ebbz
      y47=y37
      y37=y27
      y27=y17
      y17=-$a2*y27-$a3*y37-$a4*y47+$b1*z17+$b2*z27+$b3*z37+$b4*z47
      $ebbf=y17
      N7=M7
      ENDIF
#LOCAL REAL z18,z28,z38,z48
#LOCAL REAL y18,y28,y38,y48
#LOCAL REAL T8
#LOCAL INTEGER M8,N8

      T8=1.0/$fs
      M8=FLOOR(TIME/T8)
      IF(M8>N8) THEN
      z48=z38
      z38=z28
      z28=z18
      z18=$eccz
      y48=y38
      y38=y28
      y28=y18
      y18=-$a2*y28-$a3*y38-$a4*y48+$b1*z18+$b2*z28+$b3*z38+$b4*z48
      $eccf=y18
      N8=M8
      ENDIF
#LOCAL REAL z21,z22,z23
#LOCAL REAL y21,y22,y23
#LOCAL REAL T2
#LOCAL INTEGER M2,N2

      T2=1.0/$fs
      M2=FLOOR(TIME/T2)
      IF(M2>N2) THEN
      z23=z22
      z22=z21
      z21=$ecz
      y23=y22
      y22=y21
      y21=-$a2*y22-$a3*y23+$b1*z21+$b2*z22+$b3*z23
      $ecf=y21
      N2=M2
      ENDIF
#LOCAL REAL z1t,z2t,z3t
#LOCAL REAL y1t,y2t,y3t
#LOCAL REAL Tt
#LOCAL INTEGER Mt,Nt

      Tt=1.0/$fs
      Mt=FLOOR(TIME/Tt)
      IF(Mt>Nt) THEN
      z3t=z2t
      z2t=z1t
      z1t=$a1t
      y3t=y2t
      y2t=y1t
      y1t=-$a2*y2t-$a3*y3t+$b1*z1t+$b2*z2t+$b3*z3t
      $a1tt=y1t
      Nt=Mt
      ENDIF
#LOCAL REAL z1b,z2b,z3b
#LOCAL REAL y1b,y2b,y3b
#LOCAL REAL Tb
#LOCAL INTEGER Mb,Nb

      Tb=1.0/$fs
      Mb=FLOOR(TIME/Tb)
      IF(Mb>Nb) THEN
      z3b=z2b
      z2b=z1b
      z1b=$a1b
      y3b=y2b
      y2b=y1b
      y1b=-$a2*y2b-$a3*y3b+$b1*z1b+$b2*z2b+$b3*z3b
      $a1tb=y1b
      Nb=Mb
      ENDIF
#LOCAL REAL z19,z29,z39
#LOCAL REAL y19,y29,y39
#LOCAL REAL T9
#LOCAL INTEGER M9,N9

      T9=1.0/$fs
      M9=FLOOR(TIME/T9)
      IF(M9>N9) THEN
      z39=z29
      z29=z19
      z19=$a19
      y39=y29
      y29=y19
      y19=-$a2*y29-$a3*y39+$b1*z19+$b2*z29+$b3*z39
      $a1t9=y19
      N9=M9
      ENDIF
