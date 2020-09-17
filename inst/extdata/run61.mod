;; 1. Based on: 60
;; 2. Description:
;;    Implementing final frem60.dir model
;; 3. Label:
;;    DAT-1, TYPE=1
$PROBLEM    run 1
$INPUT      NO ID STUDYID TAD TIME DAY AMT RATE ODV DV EVID BLQ DOSE
            FOOD FORM TYPE WT HT AGE SEX RACE ETHNIC GENO SMOK AST ALT
            BILI BMI CRCL NCI NCIL RACEL RACEL1 RACEL2 RACEL3 GENO1
            GENO2 GENO3 GENO4 FREMTYPE
$DATA      frem60.dir/frem_dataset.dta IGNORE=@
$SUBROUTINE ADVAN2 TRANS2
$PK
;;; FRELGENO-DEFINITION START
IF(GENO.EQ.2.0000E+00) FRELGENO = 1  ; Most common
IF(GENO.EQ.3.0000E+00) FRELGENO = ( 1 + THETA(14))
IF(GENO.EQ.4.0000E+00) FRELGENO = ( 1 + THETA(15))
IF(GENO.EQ.1.0000E+00) FRELGENO = ( 1 + THETA(16))
;;; FRELGENO-DEFINITION END


;;; FRELFORM-DEFINITION START
IF(FORM.EQ.1.0000E+00) FRELFORM = 1  ; Most common
IF(FORM.EQ.0.0000E+00) FRELFORM = ( 1 + THETA(13))
;;; FRELFORM-DEFINITION END

;;; FREL-RELATION START
FRELCOVTIME = FRELFORM
FRELCOV = FRELGENO
;;; FREL-RELATION END


;;; CLFOOD-DEFINITION START
IF(FOOD.EQ.1.0000E+00) CLFOOD = 1  ; Most common
IF(FOOD.EQ.0.0000E+00) CLFOOD = ( 1 + THETA(12))
;;; CLFOOD-DEFINITION END

CLWT = (WT/75)**THETA(2)

IF(GENO.EQ.2) CLGENO = 1
IF(GENO.EQ.1) CLGENO = (1+THETA(9))
IF(GENO.EQ.3) CLGENO = (1+THETA(10))
IF(GENO.EQ.4) CLGENO = (1+THETA(11))

;;; CL-RELATION START
CLCOVTIME = CLFOOD
CLCOV=CLWT*CLGENO
;;; CL-RELATION END

VWT =  (WT/75)**THETA(3)
VCOV = VWT

TVFREL  = THETA(1)*FRELCOV
TVCL    = THETA(4)*CLCOV
TVV     = THETA(5)*VCOV
TVMAT   = THETA(6)
TVD1    = THETA(7)
TVRUV   = THETA(8)

MU_1  = LOG(TVRUV)
MU_2  = TVD1
MU_3  = LOG(TVFREL)
MU_4  = LOG(TVCL)
MU_5  = LOG(TVV)
MU_6  = LOG(TVMAT)

RUV   = EXP(MU_1               + ETA(1))
D1FR  = MU_2                   + ETA(2)
FREL  = FRELCOVTIME * EXP(MU_3 + ETA(3))
CL    = CLCOVTIME   * EXP(MU_4 + ETA(4))
V     = EXP(MU_5               + ETA(5))
MAT   = EXP(MU_6               + ETA(6))

D1    = MAT*(1-D1FR)
F1    = FREL
KA    = 1 / (MAT-D1)
S2    = V

     MU_7 = THETA(17)
     COV7 = MU_7 + ETA(7)
     MU_8 = THETA(18)
     COV8 = MU_8 + ETA(8)
     MU_9 = THETA(19)
     COV9 = MU_9 + ETA(9)
     MU_10 = THETA(20)
     COV10 = MU_10 + ETA(10)
     MU_11 = THETA(21)
     COV11 = MU_11 + ETA(11)
     MU_12 = THETA(22)
     COV12 = MU_12 + ETA(12)
$ERROR
CP    = A(2)*1000 / V
IPRED = LOG(CP + 0.00001)
Y     = IPRED + EPS(1) * EXP(ETA(1))

;;;FREM CODE BEGIN COMPACT
;;;DO NOT MODIFY
     IF (FREMTYPE.EQ.100) THEN
;       AGE  1
        Y = COV7 + EPS(2)
        IPRED = COV7
     END IF
     IF (FREMTYPE.EQ.200) THEN
;       BMI  1
        Y = COV8 + EPS(2)
        IPRED = COV8
     END IF
     IF (FREMTYPE.EQ.300) THEN
;       CRCL  1
        Y = COV9 + EPS(2)
        IPRED = COV9
     END IF
     IF (FREMTYPE.EQ.400) THEN
;       SEX  1
        Y = COV10 + EPS(2)
        IPRED = COV10
     END IF
     IF (FREMTYPE.EQ.500) THEN
;       RACEL_3  1
        Y = COV11 + EPS(2)
        IPRED = COV11
     END IF
     IF (FREMTYPE.EQ.600) THEN
;       RACEL_2  1
        Y = COV12 + EPS(2)
        IPRED = COV12
     END IF
;;;FREM CODE END COMPACT
$THETA  1 FIX ; 1. TVFREL
$THETA  0.75 FIX ; 2. WT coef on TVCL
$THETA  1 FIX ; 3. WT coef on TVV
$THETA  (0,12.9238) ; 4. TVCL
$THETA  (0,206.447) ; 5. TVV
$THETA  (0,6.6169) ; 6. TVMAT
$THETA  (0,0.63728,1) ; 7 TVD1
$THETA  1 FIX ; TVRUV
$THETA  (-1,0.733499) ; GENO=1 on CL
$THETA  (-1,-0.257239) ; GENO=3 on CL
$THETA  (-1,-0.519219) ; GENO=4 on CL
$THETA  (-1.00,0.24303,3.00) ; CLFOOD1
$THETA  (-1.00,-0.143368,3.00) ; FRELFORM1
$THETA  (-1.00,0.322536,3.00) ; FRELGENO1
 (-1.00,2.77068,3.00) ; FRELGENO2
 (-1.00,-0.408129,3.00) ; FRELGENO3
$THETA  44.2148 ; TV_AGE
 30.3011 ; TV_BMI
 118.823 ; TV_CRCL
 1.4465 ; TV_SEX
 0.0193571 ; TV_RACEL_3
 0.200529 ; TV_RACEL_2
$OMEGA  0.0821397  ; 1. IIV on RUV
$OMEGA  0.0001  FIX  ;         D1
$OMEGA  BLOCK(10)
 0.355118  ; 2. IIV on FREL
 0.107884 0.162604  ; 3. IIV on CL
 0.0632937 0.0540479 0.146033  ; 4. IIV on V
 -0.000997759 0.0275816 0.0511788 0.151893  ; 5. IIV on MAT
 1.35309 0.940607 1.13431 -0.116914 171.152  ;    BSV_AGE
 1.80361 1.53831 1.83275 0.180495 21.5876 39.8735  ;    BSV_BMI
 5.53285 4.99748 5.83668 1.00508 -119.198 79.9421 592.248  ;   BSV_CRCL
 -0.058362 -0.0226342 -0.0124707 -0.00391241 -0.248546 -0.12715 -0.107272 0.247506  ;    BSV_SEX
 -0.00114649 0.00132732 -0.0046923 0.00394207 0.0319504 -0.0576348 -0.295547 -0.00173723 0.019015  ; BSV_RACEL_3
 -0.00593722 -0.00209445 0.00427665 -0.00334753 -1.12644 -0.0188588 0.276756 0.00169173 -0.0038889 0.160554  ; BSV_RACEL_2
$SIGMA  0.0429072  ;     1. RUV
$SIGMA  1E-07  FIX  ;     EPSCOV
$ESTIMATION METHOD=IMPMAP AUTO=1 RANMETHOD=3P INTER NOABORT PRINT=1
            CITER=15 NOCOV=1 ISAMPEND=100 NITER=100 ISAMPLE=300 CTYPE=1
            CALPHA=0.001
$ESTIMATION METHOD=IMPMAP INTER EONLY=1 NITER=5 ISAMPLE=1000 NOCOV=0
$COVARIANCE PRINT=E UNCONDITIONAL


