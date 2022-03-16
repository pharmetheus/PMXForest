;; 1. Based on: frem60-2
;; 2. Description:
;;    Re-run of frm model
;; 3. Label:
;;    DAT-1, TYPE=1
$PROBLEM    run 1
$INPUT      NO ID STUDYID TAD TIME DAY AMT RATE ODV DV EVID BLQ DOSE
            FOOD FORM TYPE WT HT AGE SEX RACE ETHNIC GENO SMOK AST ALT
            BILI BMI CRCL NCI NCIL RACEL RACEL1 RACEL2 RACEL3 GENO1
            GENO2 GENO3 GENO4 FREMTYPE
$DATA      frem60-2.dir/frem_dataset.dta IGNORE=@
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
;       RACEL2  1
        Y = COV11 + EPS(2)
        IPRED = COV11
     END IF
;;;FREM CODE END COMPACT
$THETA  1 FIX ; 1. TVFREL
$THETA  0.75 FIX ; 2. WT coef on TVCL
$THETA  1 FIX ; 3. WT coef on TVV
$THETA  (0,12.9203) ; 4. TVCL
$THETA  (0,206.073) ; 5. TVV
$THETA  (0,6.62191) ; 6. TVMAT
$THETA  (0,0.637306,1) ; 7 TVD1
$THETA  1 FIX ; TVRUV
$THETA  (-1,0.75811) ; GENO=1 on CL
$THETA  (-1,-0.250755) ; GENO=3 on CL
$THETA  (-1,-0.528443) ; GENO=4 on CL
$THETA  (-1.00,0.242054,3.00) ; CLFOOD1
$THETA  (-1.00,-0.143198,3.00) ; FRELFORM1
$THETA  (-1.00,0.31564,3.00) ; FRELGENO1
 (-1.00,2.70051,3.00) ; FRELGENO2
 (-1.00,-0.394359,3.00) ; FRELGENO3
$THETA  44.2081 ; TV_AGE
 30.2923 ; TV_BMI
 118.796 ; TV_CRCL
 1.44678 ; TV_SEX
 0.200556 ; TV_RACEL2
$OMEGA  0.0844544  ; 1. IIV on RUV
$OMEGA  0.0001  FIX  ;         D1
$OMEGA  BLOCK(9)
 0.353835  ; 2. IIV on FREL
 0.108382 0.162928  ; 3. IIV on CL
 0.0640972 0.0533095 0.144797  ; 4. IIV on V
 -0.000893337 0.02743 0.0505315 0.150268  ; 5. IIV on MAT
 1.38069 0.960561 1.09519 -0.075329 171.152  ;    BSV_AGE
 1.77461 1.51175 1.81682 0.161892 21.5876 39.8734  ;    BSV_BMI
 5.51087 4.99023 5.89499 0.981828 -119.198 79.9419 592.247  ;   BSV_CRCL
 -0.0537116 -0.0183754 -0.0128324 -0.00250836 -0.248576 -0.12717 -0.107278 0.247506  ;    BSV_SEX
 -0.00566168 -0.00213343 0.00460884 -0.00451285 -1.12644 -0.0188651 0.276718 0.00169106 0.160555  ; BSV_RACEL2
$SIGMA  0.0428509  ;     1. RUV
$SIGMA  1E-07  FIX  ;     EPSCOV
$ESTIMATION METHOD=IMPMAP AUTO=1 RANMETHOD=3P INTER NOABORT PRINT=1
            CITER=15 NOCOV=1 ISAMPEND=100 NITER=50 ISAMPLE=300 CTYPE=1
            CALPHA=0.001
$ESTIMATION METHOD=IMPMAP INTER EONLY=1 NITER=5 ISAMPLE=1000 NOCOV=0
$COVARIANCE PRINT=E UNCONDITIONAL

