$PROBLEM    Time-to-event model example Weibull distribution
$INPUT      ID TIME DV EVID DOSE EXPO AGE
$DATA       tte_data1.dat
            IGNORE=@
$SUBROUTINE ADVAN=13 TOL=9
$MODEL      COMP=(HAZARD)
$PK

   LAM= THETA(1)*EXP(ETA(1))
   SHP=THETA(2)

;--------------Covariate relationships---------------
;--------------SCM compatible code-------------------
   TVRF =  (AGE/50)**THETA(3)
   TVEFF = THETA(4)*EXPO
   RFCOV = 0

   TVRF = RFCOV+TVRF
   RF   = 1*TVRF+TVEFF



$DES


   DEL= 1E-12

;----------hazard-----------------------------------

   DADT(1)=LAM*EXP(RF)*SHP*(LAM*T+DEL)**(SHP-1)

;----------TTE Model------------------------------

$ERROR

  CHZ = A(1)
  SURX = EXP(-CHZ)            ;survival probability

  DELX = 1E-12

  HAZNOW = LAM*SHP*EXP(RF)*(LAM*TIME+DELX)**(SHP-1)

  Y=SURX                      ;censored event (prob of survival)
  IF(DV.EQ.1)  Y=SURX*HAZNOW  ;prob density function of event

$THETA  (0,5E-03)    ; LAMBDA; 1
$THETA  (0.2,1.3)    ; SHAPE; 2  1 FIX for exponential distr
$THETA  (0,0.5,2)    ; AGE EFFECT
$THETA  (-1,-0.01,1) ; EXPOSURE EFFECT

$OMEGA  0  FIX
$ESTIMATION MAXEVAL=9999 METHOD=0 LIKE SIGL=9 NSIG=3 PRINT=1
$COVARIANCE MATRIX=R UNCONDITIONAL PRINT=E
$TABLE      ID TIME SURX CHZ EVID AGE EXPO DOSE NOPRINT ONEHEADER FILE=tte_weibull_sdtab

