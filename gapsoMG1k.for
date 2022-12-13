C***********************************************************************
C
C      MULTIOBJECTIVE BUFFER & SERVICE (MU) ALLOCATION IN
C      M/G/1/K QUEUEING NETWORKS + \SUM P_K + MOPSO
C      VERSION: 2020.03.28.1
C
C***********************************************************************
       PROGRAM MUBAPGAPSO_MG1K
C***********************************************************************
      INTEGER IARC, IQEU
      PARAMETER (IARC=780, IQEU=40)
C      IARC      MAXIMUM NUMBER OF ARCS (IARC=IQEU*(IQEU-1)/2)
C      IQEU      MAXIMUM NUMBER OF QUEUES
      COMMON/ONE/XLEFF,RATE(40),RH,RIC(40),RHO,XL(40),NQ,AS(40),ARH(40),
     + NARCS,NS(780),NF(780),STEST(40),SLAM(40),RP(780),ETEST(40),NRUN,
     + BETA(40),GAM(40),M,PC(40),PZ(40),EXL,TCAP,AMU,THTPUT,TOT,PENALTY
      COMMON/TWO/X(40),Y(40),
     + S(40),FY,N,KOUNTS,LIN,NDRV,DIRECT(40,40),BEFORE(40),
     + FIRST(40),SECND(40)
      COMMON/THREE/W(40),FX
      INTEGER NQ,NSNOD,SOUR(40),ENDD(40),NEND,I,J
      REAL WINT(40)
      REAL CPU_BEG,CPU_END
C      INTEGER ORDER
C      REAL BESTW(40),BESTFX,BESTTOT
      CHARACTER (LEN=100) FMT
C
C      AMU     COST PARAMETER IN THE RELAXED PROBLEM
C      ARH     SERVICE RATE AT HOLDING NODE
C      AS      SQUARED COEFICIENT OF VARIATION OF THE SERVICE TIME
C      BEFORE
C      BESTFX  BEST OBJECTIVE FUNCTION VALUE
C      BESTTOT BEST TOTAL THROUGHPUT
C      BETA
C      DIRECT
C      ENDD
C      ETEST
C      EXL     EXPECTED QUEUE LENGTH
C      FIRST
C      FMT     FORMAT SPECIFIER
C      FX
C      FY
C      GAM
C      I
C      J
C      K
C      KOUNTS
C      LIN
C      M       FINITE QUEUE BUFFER SIZE
C      N
C      NARCS
C      NDRV
C      NEND
C      NF      END NODE FOR THE ARCS
C      NQ      NUMBER OF FINITE QUEUES
C      NRUN
C      NS      START NODE FOR THE ARCS
C      NSNOD   NUMBER OF STARTING NODES
C      PC      PROBABILITY QUEUE IS AT CAPACITY
C      PENALTY
C      PZ      PROBABILITY OF 0 BUZY SERVERS AT A QUEUE
C      RATE    SERVICE RATE AT NODE I
C      RH      OVERAL SERVICE RATE AT HOLDING NODE
C      RHO     RHO OF QUEUE I
C      RIC     PURE BUFFER SIZE OF QUEUE I (EXCLUDING THOSE IN SERVICE)
C      RP
C      S
C      SECND
C      SLAM    INITIAL ARRIVAL RATE TO NODE I
C      SOUR    SOURCE NODES
C      STEST
C      TCAP
C      THTPUT
C      CPU_BEG TIME AT THE BEGINING OF RUNNING
C      CPU_END TIME AT THE END OF RUNNING
C      TOT
C      W
C      WINT
C      X
C      XL
C      XLEFF   EFFECTIVE ARRIVAL RATE (LAMBDA EFF)
C      XS      VECTOR X FOR MULTIOBJECTIVE OPTIMIZATION
C      FXS     VECTOR OF OBJETIVE FUNCTIONS
C      Y
C
C      OPTIMIZATION ALGORITHM DATA STRUCTURES
C
      INTEGER MAXPOP, MAXBAS, MAXOBJ, NGEN, NCYC, POPSIZE, NBASE, NOBJ
      PARAMETER (MAXPOP=2000, MAXBAS=200, MAXOBJ=3)
      REAL XS(MAXPOP,MAXBAS), FXS(MAXPOP,MAXOBJ)
C
C      MAXPOP  GA MAXIMUM POPULATION
C      MAXBAS  GA MAXIMUM NUMBER OF BASES
C      MAXOBJ  GA MAXIMUM NUMBER OF OBJECTIVE FUNCTIONS
C      POPSIZE GA POPULATION SIZE
C      NBASE   GA NUMBER OF BASES
C      NOBJ    GA NUMBER OF OBJECTIVE FUNCTIONS
C      XS      GA POPULATION
C      FXS     GA OBJECTIVE FUNCTIONS
C
C      TESTING RANDON NUMBER GENERATORS
C
C      DOUBLE PRECISION SEED0
C      REAL NO(10000), ETA, BET
C      SEED0=246813579.d0
C      ETA=8.0
CC
CC      NORMAL RANDON VARIABLE
CC
C      CALL NORM(SEED0,NO,10000)
C      WRITE(6,1) (NO(I),I=1,10000)
C    1 FORMAT(1(F10.7))
CC
CC      BETA RANDON VARIABLE - 1
CC
C      CALL UNI(SEED0,NO,10000)
C      DO 2 I=1,10000
C            IF (NO(I).LE.0.5) THEN
C                  WRITE(6,*) (2*NO(I))**(1/(ETA+1))
C            ELSE
C                  WRITE(6,*) (1/(2*(1-NO(I))))**(1/(ETA+1))
C            ENDIF
C    2 CONTINUE
CC
CC      BETA RANDON VARIABLE - 2
CC
C      WRITE(6,*) 'UNI:    BETA:'
C    4 CALL UNI(SEED0,NO,1)
C      IF (NO(1).LT.0.5) THEN
C            CALL UNI(SEED0,NO(1),1)
C            IF (NO(1).LE.0.5) THEN
C                  BET=(2*NO(1))**(1/(ETA+1))
C            ELSEIF (NO(1).LT.(1-1E-08)) THEN
C                  BET=(1/(2*(1-NO(1))))**(1/(ETA+1))
C            ELSE
C                  BET=(1/(2*(1E-08)))**(1/(ETA+1))
C            ENDIF
C            WRITE(6,*) NO(1), BET
C      ENDIF
C      IF (BET.LT.100) GOTO 4
C      STOP
C
C      INITIALIZE VARIABLES
C
      KOUNTS=0
      ICONVG=2
      STEP=1.0
      LIN=0
      DO 10 I=1,IQEU
            SOUR(I)=0
            STEST(I)=0
            SLAM(I)=0.0
            ENDD(I)=0.0
            ETEST(I)=0
   10 CONTINUE
C
C      OPEN INPUT FILE
C
C      OPEN(UNIT=5,NAME='SeNq3.txt',STATUS='OLD')
C      OPEN(UNIT=5,FILE='CONIN$',STATUS='OLD')
C
C      CREATE OUTPUT FILE
C
C      OPEN(UNIT=6,NAME='OUTPUT.TXT',STATUS='UNKNOWN')
C      OPEN(UNIT=6,FILE='CONOUT$',STATUS='UNKNOWN')
C
      READ (5,*)
      READ (5,*) NQ
      WRITE(6,*) 'NUMBER OF QUEUES, NQ'
      WRITE(6,*) NQ
C
      READ (5,*)
      READ (5,*) NARCS
      WRITE(6,*) 'NUMBER OF ARCS, NARCS'
      WRITE(6,*) NARCS
C
      READ (5,*)
      WRITE(6,*) 'STARTING & ENDING NODES, AND ',
     +  'ROUTING PROBABILITIES, NS(I), NF(I), RP(I)'
      DO 20 I=1,NARCS
            READ(5,*) NS(I),NF(I),RP(I)
            WRITE(6,*) NS(I),NF(I),RP(I)
   20 CONTINUE
C
C      ENTER SOURCE NODES
C
      READ (5,*)
      READ (5,*) NSNOD
      WRITE(6,*) 'NUMBER OF SOURCE NODES, NSNOD'
      WRITE(6,*) NSNOD
C
      READ (5,*)
      WRITE(6,*) 'SOURCE NODES AND ARRIVAL RATE, SOUR, SLAM(I)'
      THTPUT=0.0
      DO 30 I=1,NSNOD
C            READ SOURCE NODE AND ARRIVAL RATE
            READ (5,*) SOUR(I), SLAM(SOUR(I))
            WRITE(6,*) SOUR(I), SLAM(SOUR(I))
C            UPDATE THRESHOLD THROUGHPUT
            THTPUT=THTPUT+SLAM(SOUR(I))
C            MARK SOURCE NODES
            STEST(SOUR(I))=1
   30 CONTINUE
C
C      ENTER END NODES
C
      READ (5,*)
      READ (5,*) NEND
      WRITE(6,*) 'NUMBER OF END NODES, NEND'
      WRITE(6,*) NEND
C
      READ (5,*)
      READ (5,*) (ENDD(J),J=1,NEND)
      WRITE(6,*) 'END NODES, ENDD'
      WRITE(6,*) (ENDD(J),J=1,NEND)
      DO 40 I=1,NEND
C            MARK END NODES
            ETEST(ENDD(I))=1
   40 CONTINUE
C
C      ENTER SERVICE RATES, SQUARED COEFF VARIATION
C      AND TOTAL CAPACITY K
C
      READ (5,*)
      WRITE(6,*) 'SERVICE RATES, COEFF VARIATION,',
     + 'AND TOT CAPACITY, RATE(I), AS(I), W(I)'
      DO 50 I=1,NQ
            READ (5,*) RATE(I),AS(I),W(I)
            ARH(I)=RATE(I)
            WRITE(6,*) RATE(I),AS(I),W(I)
   50 CONTINUE
C
C      EVALUATE SUM OF BLOCKING PROBABILITIES SUM PK
C
      DO 60 I=1,NQ
            WINT(I)=W(I)
            WINT(NQ+I)=RATE(I)
   60 CONTINUE
      WRITE(6,*)
      WRITE(6,*) 'CALLING FUN(WINT,FX):'
      CALL FUN(WINT,FX,.FALSE.)
      FXS(1,1)=0.0
      FXS(1,2)=0.0
      DO 62 I=1,NQ
            FXS(1,1)=FXS(1,1)+WINT(I)
            FXS(1,2)=FXS(1,2)+WINT(NQ+I)
   62 CONTINUE
      FXS(1,3)=FX
C
C      PRINT RESULTS
C
      WRITE(6,*) 'RESULTS:'
      WRITE(6,*) 'SUM_K    SUM_MU     SUM_PK'
      WRITE(6,80) INT(FXS(1,1)), FXS(1,2), FXS(1,3)
   80 FORMAT(1(I6),(F10.4),(F10.6))
C
C      STOP HERE
C
C      STOP
C
C      DON'T GENERATE 3-D FIGURE
C
      GOTO 160
C
C      GENERATE 3-D FIGURE
C
      WRITE(6,*)
      WRITE(6,*) '3-D FIGURE (K VS MU VS THETA):'
      DO 150 I=1, 20
C            SET SERVICE RATE
            RATE(1)=FLOAT(I)/2.+1E-06
            DO 140 J=1,40
                  WINT(1)=J
                  WINT(2)=RATE(1)
C                  WRITE(6,*)
C                  WRITE(6,*) 'CALLING FUN(WINT,FX):'
                  CALL FUN(WINT,FX,.FALSE.)
C                  WRITE(6,*) 'RESULTS:'
C                  WRITE(6,*) 'SUM_K, LAMBDA, THETA
                  WRITE(6,130) INT(WINT(1)), WINT(2), TOT
  130             FORMAT(1(I4),(F10.4),(F10.6))
  140       CONTINUE
  150 CONTINUE
      STOP
  160 CONTINUE
C
C      INITIALIZE PARAMETERS AND OPTIMIZE BY NSGA-II
C
C      NUMBER OF GENERATIONS
C
      NGEN=4000
C
C      POPULATION SIZE
C
      POPSIZE=400
C
C      USE NSGA-II
C
      NBASE=2*NQ
      NOBJ=3
      CALL CPU_TIME(CPU_BEG)
      WRITE(6,*)
      WRITE(6,*) 'CALLING MOGA:'
      CALL NSGAII(NGEN,POPSIZE,NBASE,NOBJ,XS,FXS)
      CALL CPU_TIME(CPU_END)
      WRITE(6,*) 'NSGAII RESULT AFTER'
      WRITE(6,*) 'GENERATIONS', NGEN
      WRITE(6,*) 'POPULATION', POPSIZE
      WRITE(6,*) 'CPU(S)', CPU_END-CPU_BEG
      WRITE(6,*) 'XS FXS'
      DO 180 I=1, POPSIZE
            WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3))")') NBASE, NOBJ
            WRITE(6,FMT) ((XS(I,J)),J=1,NBASE), (FXS(I,J),J=1,NOBJ)
  180 CONTINUE
C
C      INITIALIZE PARAMETERS AND OPTIMIZE BY MOPSO
C
C      NUMBER OF CYCLES
C
      NCYC=4000
C
C      POPULATION SIZE
C
      POPSIZE=400
C
      NBASE=2*NQ
      NOBJ=3
      CALL CPU_TIME(CPU_BEG)
      WRITE(6,*)
      WRITE(6,*) 'CALLING MOPSO:'
      CALL MOPSO(NCYC,POPSIZE,NBASE,NOBJ,XS,FXS)
      CALL CPU_TIME(CPU_END)
      WRITE(6,*) 'MOPSO RESULT AFTER'
      WRITE(6,*) 'CYCLES', NCYC
      WRITE(6,*) 'POPULATION', POPSIZE
      WRITE(6,*) 'CPU(S)', CPU_END-CPU_BEG
      WRITE(6,*) 'XS FXS'
      DO 190 I=1, POPSIZE
            WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3))")') NBASE, NOBJ
            WRITE(6,FMT) ((XS(I,J)),J=1,NBASE), (FXS(I,J),J=1,NOBJ)
  190 CONTINUE
C
      CLOSE(UNIT=5)
      CLOSE(UNIT=6)
      END
C***********************************************************************
C
C     MULTI-OBJECTIVE PARTICLE SWARM OPTIMIZATION ALGORITHM
C
C***********************************************************************
      SUBROUTINE MOPSO(NCYC,POPSIZE,NBASE,NOBJ,XS,FXS)
C***********************************************************************
      DOUBLE PRECISION SEED0
      EXTERNAL FUN
      INTEGER MAXPOP, MAXBAS, MAXOBJ, NCYC, POPSIZE, NBASE, NOBJ
      PARAMETER (MAXPOP=2000, MAXBAS=200, MAXOBJ=3)
      INTEGER I, J, K, L, M, QSIZE
      INTEGER FSIZE(2*MAXPOP), NP(2*MAXPOP), RANK(2*MAXPOP)
      INTEGER Q(2*MAXPOP), SP(2*MAXPOP,2*MAXPOP)
      INTEGER SPSIZE(2*MAXPOP), F(2*MAXPOP,2*MAXPOP), ISORT(2*MAXPOP)
      INTEGER IMIN, IAUX, PAR1, BETTER, WORSE
      REAL MAXXS, FX, UN, NO(2)
      REAL MAXFIT, XS(MAXPOP,MAXBAS), FXS(MAXPOP,MAXOBJ)
      REAL X(MAXBAS), PT(MAXPOP,MAXBAS), QT(MAXPOP,MAXBAS)
      REAL RT(2*MAXPOP,MAXBAS), IOBJ(2*MAXPOP,MAXOBJ), IDIST(2*MAXPOP)
      REAL VEL(MAXPOP,MAXBAS),PBEST(MAXPOP,MAXBAS),GBEST(MAXPOP,MAXBAS)
      REAL WW
      REAL ETA, QLIM, DL(MAXPOP), MAXID, SIGL, SUMM, SUM2, DELTAL
      CHARACTER (LEN=20) FILENAME
      INTEGER IDX, LDL, NDOM
      LOGICAL IPRINT, SNAP, CONVDEL, CONVPS, CONVDOM
      PARAMETER (ETA=8.0, LDL=40, DELTAL=0.02)
      CHARACTER (LEN=100) FMT
C      INERTIA WEIGHT
      PARAMETER (WW=0.4)
C
C     TO PRINT RESULTS IPRINT=.TRUE., OTHERWISE, IPRINT=.FALSE.
C
      IPRINT=.FALSE.
C
C     F       SETS FI
C     FSIZE   SET FI SIZE
C     FX      THROUGHPUT
C     MAXXS   UPPER BOUND FOR XS VECTOR
C     IS      SET IS
C     NCYC    NUMBER OF CYCLES
C     POPSIZE POPULATION SIZE
C     MAXPOP  MAXIMUM POPULATION SIZE
C     NBASE   NUMBER OF BASES IN THE GENETIC CODE
C     NOBJ    NUMBER OF OBJECTIVE FUNCTIONS
C     MAXBAS  MAXIMUM NUMBER OF BASES
C     MAXOBJ  MAXIMUM NUMBER OF OBJECTIVES
C     PT      SET PT
C     QT      SET QT
C     Q       SET Q
C     QSIZE   SET Q SIZE
C     RANK    SET OF RANKS
C     NP      SET NP
C     SP      SET SP
C     SPSIZE  SET SP SIZE
C     RT      SET RT
C     IOBJ    OBJECTIVE FUNCTIONS OF SET RT
C     ISORT   SET SORTED I
C     IAUX    AUXILIARY VARIABLE FOR SORTING
C     IMIN    AUXILIARY VARIABLE FOR SORTING
C     SEED0   SEED FOR RANDOM NUMBER GENERATION
C     FUN     IMPLEMENTS THE GENERALIZED EXPANSION METHOD
C     XS      VECTOR OF CAPACITIES
C     FXS     VECTOR OF FUNCTION VALUES
C     UN      UNIFORMLY DISTRIBUTED RANDOM NUMBER
C     CROSTYP GA TYPE OF CROSSOVER (0, NO; 1, UNIFORM; 2, FITNESS ORIENTED)
C     RATEMUT GA MUTATION RATE (WORKS BETTER BETWEEN 1%-2%)
C     PBEST   POPULATION BEST POSITION
C     GBEST   POPULATION GLOBAL LEADER
C     WW      INERTIA WEIGHT
C     Q       SET Q
C
      IF (IPRINT) WRITE(6,*) 'MOPSO'
C
C     SET SEED FOR RANDOM NUMBER GENERATION
C
      SEED0=246813579.d0
C     SEED0=135792468.d0
C
C     GENERATE INITIAL POPULATION AND SPEED
C
      MAXXS=60
C      UN=1.0
      DO 20 I=1, POPSIZE
C          CALL UNI(SEED0,NO,1)
C          UN=NO(1)
          DO 10 J=1, NBASE/2
C
C             PT'S
C
C              CALL UNI(SEED0,NO,1)
C              UN=NO(1)
C              PT(I,J)=CEILING(UN*MAXXS)
              PT(I,J)=XS(I,J)
C             WRITE(6,*) 'CEILING(UN*MAXXS)', CEILING(UN*MAXXS)
C              CALL UNI(SEED0,NO,1)
C              UN=NO(1)
C              PT(I,NBASE/2+J)=(UN*MAXXS)
              PT(I,NBASE/2+J)=XS(I,NBASE/2+J)
C             WRITE(6,*) '(UN*MAXXS)', (UN*MAXXS)
C
C             QT'S
C
C              CALL UNI(SEED0,NO,1)
C              UN=NO(1)
C              QT(I,J)=CEILING(UN*MAXXS)
              QT(I,J)=XS(I,J)
C             WRITE(6,*) 'CEILING(UN*MAXXS)', CEILING(UN*MAXXS)
C              CALL UNI(SEED0,NO,1)
C              UN=NO(1)
C              QT(I,NBASE/2+J)=(UN*MAXXS)
              QT(I,NBASE/2+J)=XS(I,NBASE/2+J)
C             WRITE(6,*) '(UN*MAXXS)', (UN*MAXXS)
C
C             SPEED
C
              VEL(I,J)=0.0
              VEL(I,NBASE/2+J)=0.0
C
   10     CONTINUE
          IF (IPRINT) THEN
              WRITE(FMT,'("(A,I0,A,",I0,"(F10.4))")') NBASE
              WRITE(6,FMT) ' PT(',I,'): ',(PT(I,J),J=1,NBASE)
              WRITE(FMT,'("(A,I0,A,",I0,"(F10.4))")') NBASE
              WRITE(6,FMT) ' QT(',I,'): ',(QT(I,J),J=1,NBASE)
              WRITE(FMT,'("(A,I0,A,",I0,"(F10.4))")') NBASE
              WRITE(6,FMT) ' VEL(',I,'): ',(VEL(I,J),J=1,NBASE)
          ENDIF
   20 CONTINUE
C
C     GENERAL INITIALIZATIONS
C
      IDX=47
      CONVDEL=.TRUE.
      CONVPS=.TRUE.
      CONVDOM=.TRUE.
      DO 30 I=1, LDL
          DL(I)=0.0
   30 CONTINUE
C
C     ITERATE FOR NCYC CYCLES
C
      IF (IPRINT) WRITE(6,*) "      CYC      BEST_SOLUTION"
      I=0
   40 CONTINUE
C         IF (I.GE.2176.AND.I.LE.2180) THEN
C             IPRINT=1
C         ELSE
C             IPRINT=0
C         ENDIF
          IF (IPRINT) WRITE(6,*) 'MOPSO: GEN', I
C
C         DEFINE SNAPSHOTS
C
          SNAP=(I.EQ.0.OR.I.EQ.10.OR.I.EQ.100.OR.I.EQ.NCYC)
C
C         UPDATE BEST POSITIONS OF PARTICLES PBEST (EQUAL TO SELECTED POPULATION PT)
C
          IF (IPRINT) WRITE(6,*) 'PBEST    '
          DO 52 J=1, POPSIZE
            DO 50 K=1, NBASE
                PBEST(J,K)=PT(J,K)
   50       CONTINUE
            IF (IPRINT) THEN
                WRITE(FMT,'("(A,I0,A,",I0,"(E10.3))")') NBASE
                WRITE(6,FMT) ' PBEST(',J,'): ',(PBEST(J,K),K=1,NBASE)
            ENDIF
   52     CONTINUE
C
C         COMBINE POPULATION(I) AND POPULATION(I+1)
C
          DO 70 J=1, POPSIZE
              DO 60 K=1, NBASE
C                 RT(2*J-1,K)=PT(J,K)
C                 RT(2*J,K)=QT(J,K)
                  RT(J,K)=PT(J,K)
                  RT(J+POPSIZE,K)=QT(J,K)
   60         CONTINUE
   70     CONTINUE
C
C         COMPUTE OBJECTIVE FUNCTIONS
C
          IF (IPRINT) WRITE(6,*) 'RT  IOBJ'
          DO 90 J=1,2*POPSIZE
              DO 80 K=1, NBASE/2
                  X(K)=RT(J,K)
                  X(NBASE/2+K)=RT(J,NBASE/2+K)
   80         CONTINUE
C             WRITE(6,*) 'MOPSO: CALLING FUN(RT,FX):'
              CALL FUN(X,FX,IPRINT)
C             UPDATE INPUT VECTOR
              DO 82 K=1, NBASE/2
                  RT(J,K)=X(K)
                  RT(J,NBASE/2+K)=X(NBASE/2+K)
   82         CONTINUE
C             COMPUTING OBJECTIVE FUNCTIONS
              IOBJ(J,1)=0.0
              IOBJ(J,2)=0.0
              DO 84 K=1, NBASE/2
                  IOBJ(J,1)=IOBJ(J,1)+X(K)
                  IOBJ(J,2)=IOBJ(J,2)+X(NBASE/2+K)
   84         CONTINUE
              IOBJ(J,3)=FX
              IF (IPRINT) THEN
               WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3))")') NBASE,NOBJ
               WRITE(6,FMT) (RT(J,K),K=1,NBASE),(IOBJ(J,K),K=1,NOBJ)
              ENDIF
   90     CONTINUE
C
C         PERFORM FAST NONDOMINATED SORTING
C
          DO 100 J=1, 2*POPSIZE
              FSIZE(J)=0
  100     CONTINUE
          DO 140 J=1, 2*POPSIZE
              SPSIZE(J)=0
              NP(J)=0
              DO 130 K=1, 2*POPSIZE
C                 DOES J DOMINATE K?
                  BETTER=0
                  WORSE=0
                  DO 110 L=1, NOBJ
C                     WRITE(6,*) 'IOBJ(J), IOBJ(K)', IOBJ(J,L),IOBJ(K,L)
                      IF (IOBJ(J,L).LT.IOBJ(K,L)) THEN
                          BETTER=BETTER+1
                      ELSE IF (IOBJ(J,L).GT.IOBJ(K,L)) THEN
                          WORSE=WORSE+1
                      ENDIF
  110             CONTINUE
C                 WRITE(6,*) 'BETTER', BETTER, ' WORSE', WORSE
                  IF (BETTER.GT.0.AND.WORSE.EQ.0) THEN
C                     WRITE(6,*) J,' DOMINATES', K, '->UPDATE SP'
                      SPSIZE(J)=SPSIZE(J)+1
                      SP(J,SPSIZE(J))=K
                  ELSE
C                     DOES K DOMINATE J?
                      BETTER=0
                      WORSE=0
                      DO 120 L=1, NOBJ
C                         WRITE(6,*) 'IOBJ(K), IOBJ(J)', IOBJ(K,L),IOBJ(J,L)
                          IF (IOBJ(K,L).LT.IOBJ(J,L)) THEN
                              BETTER=BETTER+1
                          ELSE IF (IOBJ(K,L).GT.IOBJ(J,L)) THEN
                              WORSE=WORSE+1
                          ENDIF
  120                 CONTINUE
C                     WRITE(6,*) 'BETTER', BETTER, ' WORSE', WORSE
                      IF ((BETTER.GT.0).AND.(WORSE.EQ.0)) THEN
C                         WRITE(6,*) K,' DOMINATES', J, '->UPDATE NP'
                          NP(J)=NP(J)+1
                      ENDIF
                  ENDIF
  130         CONTINUE
              IF (IPRINT) THEN
                  IF (SPSIZE(J).GT.0) THEN
                   WRITE(FMT,'("(A,I0,A,",I0,"I4)")') SPSIZE(J)
                   WRITE(6,FMT)' SP(',J, '): ',((SP(J,K)),K=1,SPSIZE(J))
                  ENDIF
              ENDIF
              IF (NP(J).EQ.0) THEN
                  RANK(J)=1
                  FSIZE(1)=FSIZE(1)+1
                  F(1,FSIZE(1))=J
              ENDIF
  140     CONTINUE
C
          IF (IPRINT) THEN
              WRITE(6,*) 'NP: '
              WRITE(6,*) (NP(J),J=1,2*POPSIZE)
              WRITE(FMT,'("(A,",I0,"I4)")') FSIZE(1)
              WRITE(6,FMT) ' F(1): ', (F(1,K),K=1,FSIZE(1))
          ENDIF
C
          J=1
  150     CONTINUE
              QSIZE=0
C             WRITE(6,*)'P',(F(J,K),K=1,FSIZE(J))
              DO 160 K=1, FSIZE(J)
C                 P=F(J,K)
C                 WRITE(6,*)'Q',(SP(F(J,K),L),L=1,SPSIZE(F(J,K)))
                  DO 154 L=1, SPSIZE(F(J,K))
                      NP(SP(F(J,K),L))=NP(SP(F(J,K),L))-1
                      IF (NP(SP(F(J,K),L)).EQ.0) THEN
                          RANK(SP(F(J,K),L))=J+1
                          QSIZE=QSIZE+1
                          Q(QSIZE)=SP(F(J,K),L)
                      ENDIF
  154             CONTINUE
C                 WRITE(6,*) 'NP', (NP(L),L=1,2*POPSIZE)
  160         CONTINUE
              J=J+1
              FSIZE(J)=0
              DO 170 K=1, QSIZE
                  FSIZE(J)=FSIZE(J)+1
                  F(J,FSIZE(J))=Q(K)
  170         CONTINUE
              IF (IPRINT) THEN
                  IF (FSIZE(J).GT.0) THEN
                    WRITE(FMT,'("(A,I0,A,",I0,"I4)")') FSIZE(J)
                    WRITE(6,FMT) ' F(',J,'): ', (F(J,K),K=1,FSIZE(J))
                  ENDIF
              ENDIF
          IF ((J.LT.(2*POPSIZE)).AND.(FSIZE(J).NE.0)) GOTO 150
C
C         COMPUTE CROWDING DISTANCES
C
          IF (IPRINT) WRITE(6,*) 'IDIST:'
          DO 180 J=1, 2*POPSIZE
              IDIST(J)=0.0
  180     CONTINUE
          DO 230 J=1, NOBJ
C             SORT BY J-TH OBJECTIVE
C             WRITE(6,*) 'OBJ', J
              DO 190 K=1, 2*POPSIZE
                  ISORT(K)=K
  190         CONTINUE
C             WRITE(6,*) 'IOBJ', (IOBJ(ISORT(L),J),L=1,2*POPSIZE)
              DO 210 K=1, 2*POPSIZE
                  IMIN=K
                  DO 200 L=K+1, 2*POPSIZE
                      IF (IOBJ(ISORT(L),J).LT.IOBJ(ISORT(IMIN),J)) THEN
                          IMIN=L
C                     ELSE
C                         IF (IOBJ(ISORT(L),J).EQ.IOBJ(ISORT(IMIN),J)
C     +                       .AND.IDIST(ISORT(L))
C     +                       .GT.IDIST(ISORT(IMIN))) THEN
C                             IMIN=L
C                         ENDIF
                      ENDIF
  200             CONTINUE
                  IAUX=ISORT(K)
                  ISORT(K)=ISORT(IMIN)
                  ISORT(IMIN)=IAUX
  210         CONTINUE
C
C             PRINT SNAPSHOT
C
C             IF (SNAP) THEN
C                 WRITE(6,*) 'MOPSO: SORTED IOBJ AFTER GENERATION', I
C                 WRITE(6,*) (IOBJ(ISORT(K),J),K=1,2*POPSIZE)
C             ENDIF
C             UPDATE CROWDING DISTANCES
              IDIST(ISORT(1))=1E+38
              IDIST(ISORT(2*POPSIZE))=1E+38
C             WRITE(6,*) 'FMAX', IOBJ(ISORT(1),J)
C             WRITE(6,*) 'FMIN', IOBJ(ISORT(2*POPSIZE),J)
              DO 220 K=2, (2*POPSIZE)-1
                  IDIST(ISORT(K))=IDIST(ISORT(K))+
     +                (IOBJ(ISORT(K+1),J)-IOBJ(ISORT(K-1),J))/
     +                (IOBJ(ISORT(2*POPSIZE),J)-IOBJ(ISORT(1),J))
  220         CONTINUE
C
C             PRINT SNAPSHOT
C
C             IF (SNAP) THEN
C                 WRITE(6,*) 'MOPSO: IDIST AFTER GENERATION', I
C                 WRITE(6,*) (IDIST(ISORT(K)), K=1,2*POPSIZE)
C             ENDIF
  230     CONTINUE
          IF (IPRINT) THEN
              WRITE(FMT,'("(",I0,"(E10.3))")') 2*POPSIZE
              WRITE(6,FMT) (IDIST(ISORT(K)), K=1,2*POPSIZE)
C              WRITE(6,232) (IDIST(ISORT(K)), K=1,2*POPSIZE)
C  232         FORMAT(<2*POPSIZE>(E10.3))
          ENDIF
C
C         PRINT SNAPSHOT OF POPULATION AT EACH LEVEL
C
          IF (SNAP) THEN
            IDX=IDX+1
            FILENAME='SNAP-'//CHAR(IDX)//'.TXT'
            OPEN(UNIT=7, FILE=FILENAME, STATUS='UNKNOWN')
C            WRITE(6,*) 'SNAPSHOT WRITTEN IN FILE ', FILENAME
            WRITE(7,*)
     +            'MOPSO: POPULATION AT EACH LEVEL AFTER GENERATION', I
            J=1
  240       CONTINUE
              DO 250 L=1, FSIZE(J)
                WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3),I4,E10.3)")')
     +                NBASE, NOBJ
                WRITE(7,FMT) ((RT(F(J,L),K)),K=1,NBASE),
     +            (IOBJ(F(J,L),K),K=1,NOBJ), J, IDIST(F(J,L))
  250         CONTINUE
              J=J+1
            IF (FSIZE(J).NE.0) GOTO 240
            CLOSE(UNIT=7)
          ENDIF
C
C         TRACK CONVERGENCE
C
C         COMPUTE MAXIMAL CROWDING DISTANCE
C
          MAXID=0.0
          DO 260 J=1, 2*POPSIZE
              IF (IDIST(J).GT.MAXID.AND.IDIST(J).LT.1E+38) THEN
                  MAXID=IDIST(J)
              ENDIF
  260     CONTINUE
C         WRITE(6,*) 'MAXID:', MAXID
C
C         STORE VALUE
C
          DL(MOD(I,LDL)+1)=MAXID
C         WRITE(6,*) 'POSITION:', MOD(I,LDL)+1
          IF (IPRINT) WRITE(6,*) 'DL:', (DL(J),J=1,LDL)
C
C         COMPUTE SIGL
C
          SUMM=0.0
          SUM2=0.0
          IF (I.GE.LDL) THEN
              DO 280 J=1, LDL
                  SUMM=SUMM+DL(J)
                  SUM2=SUM2+DL(J)**2
  280         CONTINUE
              SIGL=SQRT(ABS(SUM2-SUMM*SUMM/LDL)/LDL)
          ELSEIF (I.NE.0) THEN
              DO 270 J=1, I+1
                  SUMM=SUMM+DL(J)
                  SUM2=SUM2+DL(J)**2
  270         CONTINUE
              SIGL=SQRT(ABS(SUM2-SUMM*SUMM/(I+1))/(I+1))
          ELSE
              SIGL=0.0
          ENDIF
          IF (IPRINT) WRITE(6,*) 'SUM2, SUMM, SIGL:', SUM2, SUMM, SIGL
          IF (CONVDEL.AND.SIGL.LT.DELTAL.AND.I.GE.LDL) THEN
            CONVDEL=.FALSE.
C           WRITE(6,*)'MOPSO: DELTA_L BELOW THRESHOLD AT GENERATION', I
          ENDIF
C
C         USE SELECTION, AND UPDATE SPEED AND POSITIONS TO CREATE NEW POPULATION
C
C         SELECTION
C
          J=1
          K=0
          IF (IPRINT) WRITE(6,*) 'PT_INI:'
  290     IF (K+(FSIZE(J)).GT.(POPSIZE)) GOTO 320
              IF (IPRINT) WRITE(6,*) ' FSIZE(',J,')=', FSIZE(J)
              DO 310 L=1, FSIZE(J)
                  K=K+1
                  DO 300 M=1, NBASE
                      PT(K,M)=RT(F(J,L),M)
  300             CONTINUE
                  IF (IPRINT) THEN
                      WRITE(FMT,'("(",I0,"(F10.4))")') NBASE
                      WRITE(6,FMT) (PT(K,M),M=1,NBASE)
                  ENDIF
  310         CONTINUE
              J=J+1
              GOTO 290
  320     CONTINUE
          IF (IPRINT) WRITE(6,*) 'PT    '
  330     IF (K.GE.POPSIZE) GOTO 360
C             WRITE(6,*) ' FSIZE(',J,')=', FSIZE(J)
C
C             FIND LARGEST CROWDING DISTANCE
C
C             ERROR DETECTED BY MEANS OF TOM'S EXAMPLE
C             IMAX=F(J,1)
C
              IMAX=1
              DO 340 L=2, FSIZE(J)
                  IF (IDIST(F(J,IMAX)).LT.IDIST(F(J,L))) THEN
                      IMAX=L
                  ENDIF
  340         CONTINUE
C             WRITE(6,*) 'IMAX', F(J,IMAX)
C
C             DELETE THIS CROWDING DISTANCE
C
              IDIST(F(J,IMAX))=-1.0
C
C             INCLUDE ELEMENT WITH LARGEST CROWDING DISTANCE
C
              K=K+1
              DO 350 M=1, NBASE
                  PT(K,M)=RT(F(J,IMAX),M)
  350         CONTINUE
              IF (IPRINT) THEN
                  WRITE(FMT,'("(",I0,"(F10.4))")') NBASE
                  WRITE(6,FMT) (PT(K,M),M=1,NBASE)
              ENDIF
              GOTO 330
  360     CONTINUE
C
C         FILL IT UP IF NECESSARY
C
          IF (IPRINT) WRITE(6,*) 'PT_COM:'
  370     IF (K.GE.POPSIZE) GOTO 390
              K=K+1
              CALL UNI(SEED0,NO,1)
              UN=NO(1)
              PAR1=CEILING(UN*2*POPSIZE)
              DO 380 M=1, NBASE
                  PT(K,M)=RT(PAR1,M)
  380         CONTINUE
              IF (IPRINT) THEN
                  WRITE(FMT,'("(",I0,"(F10.4))")') NBASE
                  WRITE(6,FMT) (PT(K,M),M=1,NBASE)
              ENDIF
              GOTO 370
  390     CONTINUE
C
C         COMPUTE GLOBAL LEDERS GBEST (RANDOMLY SELECTED FROM BEST FRONT F(1))
C
          IF (IPRINT) WRITE(6,*) 'GBEST    '
          DO 430 J=1, POPSIZE
            CALL UNI(SEED0,NO,1)
            L=CEILING(NO(1)*FSIZE(1))
            DO 420 K=1, NBASE
                GBEST(J,K)=RT(F(1,L),K)
  420       CONTINUE
            IF (IPRINT) THEN
                WRITE(FMT,'("(A,I0,A,",I0,"(F10.4))")') NBASE
                WRITE(6,FMT) ' GBEST(',J,'): ',(GBEST(J,K),K=1,NBASE)
            ENDIF
  430     CONTINUE
C
C         COMPUTE OBJECTIVE FUNCTIONS (REDUNDANT)
C
          IF (IPRINT) WRITE(6,*) 'PT IOBJ'
          MAXFIT=0.0
          DO 470 J=1,POPSIZE
              DO 440 K=1, NBASE/2
                  X(K)=PT(J,K)
                  X(NBASE/2+K)=PT(J,NBASE/2+K)
  440         CONTINUE
C             WRITE(6,*)
C             WRITE(6,*) 'MOPSO: CALLING FUN(PT,FX):'
              CALL FUN(X,FX,IPRINT)
C             UPDATE INPUT VECTOR
              DO 450 K=1, NBASE/2
                  PT(J,K)=X(K)
                  PT(J,NBASE/2+K)=X(NBASE/2+K)
  450         CONTINUE
C             COMPUTING OBJECTIVE FUNCTIONS
              IOBJ(J,1)=0.0
              IOBJ(J,2)=0.0
              DO 460 K=1, NBASE/2
                  IOBJ(J,1)=IOBJ(J,1)+X(K)
                  IOBJ(J,2)=IOBJ(J,2)+X(NBASE/2+K)
  460         CONTINUE
              IOBJ(J,3)=FX
              IF (-IOBJ(J,3).GT.MAXFIT) THEN
                  MAXFIT=-IOBJ(J,3)
              ENDIF
              IF (IPRINT) THEN
                  WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3))")')
     +                      NBASE, NOBJ
                  WRITE(6,FMT) (PT(J,K),K=1,NBASE),(IOBJ(J,K),K=1,NOBJ)
              ENDIF
C             WRITE(6,406) (PT(J,K),K=1,NBASE),(IOBJ(J,K),K=1,NOBJ)
  470     CONTINUE
C
C         TRACK CONVERGENCE
C
C         COMPUTE NUMBER OF DOMINATED SOLUTIONS (AT PARETO SET ONLY)
C
          NDOM=0
          DO 500 K=1, MIN(POPSIZE,FSIZE(1))
C             DOES FXS DOMINATES K?
              BETTER=0
              WORSE=0
              DO 480 L=1, NOBJ
C                 WRITE(6,*) 'FXS, IOBJ(K)', FXS(1,L),IOBJ(K,L)
                  IF (FXS(1,L).LT.IOBJ(K,L)) THEN
                      BETTER=BETTER+1
                  ELSE IF (FXS(1,L).GT.IOBJ(K,L)) THEN
                      WORSE=WORSE+1
                  ENDIF
  480         CONTINUE
C             WRITE(6,*) 'BETTER', BETTER, ' WORSE', WORSE
              IF (BETTER.GT.0.AND.WORSE.EQ.0) THEN
C                 WRITE(6,*) 'FXS DOMINATES', K, '->UPDATE NDOM'
                  NDOM=NDOM+1
              ELSE
C                 DOES K DOMINATE FXS?
                  BETTER=0
                  WORSE=0
                  DO 490 L=1, NOBJ
C                     WRITE(6,*) 'IOBJ(K), FXS', IOBJ(K,L),FXS(1,L)
                      IF (IOBJ(K,L).LT.FXS(1,L)) THEN
                          BETTER=BETTER+1
                      ELSE IF (IOBJ(K,L).GT.FXS(1,L)) THEN
                          WORSE=WORSE+1
                      ENDIF
  490             CONTINUE
C                 WRITE(6,*) 'BETTER', BETTER, ' WORSE', WORSE
                  IF (CONVDOM.AND.BETTER.GT.0.AND.WORSE.EQ.0) THEN
                      CONVDOM=.FALSE.
                      WRITE(6,*)
     +                     'MOPSO: WARNING: FXS HAS BEEN DOMINATED'
                  ENDIF
              ENDIF
  500     CONTINUE
C
C         PRINT OUT CONVERGENCE RESULTS
C
C         WRITE(6,*) I, SIGL, NDOM
C
C         WARN WHEN CONVERGENCE IS FIRST REACHED
C
          IF (CONVPS.AND.NDOM.EQ.0) THEN
            CONVPS=.FALSE.
C           WRITE(6,*)'MOPSO: TARGET IN PARETO SET AT GENERATION', I
          ENDIF
C
C         PRINT SNAPSHOT OF THE RESULT (PARETO-SET ONLY) AFTER GENERATION I
C
          IF (SNAP) THEN
              IDX=IDX+1
              FILENAME='SNAP-'//CHAR(IDX)//'.TXT'
              OPEN(UNIT=7, FILE=FILENAME, STATUS='UNKNOWN')
C              WRITE(6,*) 'SNAPSHOT WRITTEN IN FILE ', FILENAME
              WRITE(7,*) 'MOPSO: RESULT AFTER GENERATION', I
              DO 510 J=1, MIN(POPSIZE,FSIZE(1))
                WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3))")')NBASE,NOBJ
                WRITE(7,FMT) ((PT(J,K)),K=1,NBASE),(IOBJ(J,K),K=1,NOBJ)
  510         CONTINUE
              CLOSE(UNIT=7)
          ENDIF
C
C         UPDATE SPEED AND POSITIONS TO CREATE NEW POPULATION
C
          IF (IPRINT) WRITE(6,*) 'QT'
          DO 530 J=1, POPSIZE
              DO 520 K=1, NBASE
                  CALL UNI(SEED0,NO,2)
                  VEL(J,K)=WW*VEL(J,K)+
     +               NO(1)*(PBEST(J,K)-PT(J,K))+
     +               NO(2)*(GBEST(J,K)-PT(J,K))
                  QT(J,K)=PT(J,K)+VEL(J,K)
C
C                 ADJUSTMENT TO INTEGRALITY, NON-NEGATIVITY, AND OVERFLOW
C
                  IF (K.LE.NBASE/2) THEN
                     QT(J,K)=CEILING(QT(J,K))
                     QLIM=100.0
                     IF (QT(J,K).GT.QLIM) THEN
                         QT(J,K)=QLIM-ABS(QT(J,K)-QLIM)
C                         WRITE(6,*) K, '-TH CORRECTED TO:', QT(J,K)
                     ENDIF
                     QLIM=1.0
                     IF (QT(J,K).LT.QLIM) THEN
                         QT(J,K)=QLIM+ABS(QT(J,K)-QLIM)
C                         WRITE(6,*) K, '-TH CORRECTED TO:', QT(J,K)
                     ENDIF
                  ELSE
                     QLIM=50.0
                     IF (QT(J,K).GT.QLIM) THEN
                         QT(J,K)=QLIM-ABS(QT(J,K)-QLIM)
C                         WRITE(6,*) K, '-TH CORRECTED TO:', QT(J,K)
                     ENDIF
                     QLIM=0.0
                     IF (QT(J,K).LT.QLIM) THEN
                         QT(J,K)=QLIM+ABS(QT(J,K)-QLIM)
C                         WRITE(6,*) K, '-TH CORRECTED TO:', QT(J,K)
                     ENDIF
                  ENDIF
  520         CONTINUE
              IF (IPRINT) THEN
                 WRITE(FMT,'("(A,I0,A,",I0,"(F10.4))")') NBASE
                 WRITE(6,FMT) ' QT(',J,'): ',(QT(J,K),K=1,NBASE)
                 WRITE(6,FMT) ' VEL(',J,'): ',(VEL(J,K),K=1,NBASE)
              ENDIF
  530     CONTINUE
          I=I+1
      IF (I.LE.NCYC) GOTO 40
C
C     FINAL RESULT
C
      IF (IPRINT) THEN
          WRITE(6,*) 'MOPSO RESULT AFTER'
          WRITE(6,*) 'CYCLES', NCYC
          WRITE(6,*) 'POPULATION', POPSIZE
          WRITE(6,*) 'XS FXS'
      ENDIF
      DO 560 I=1, POPSIZE
          DO 540 J=1, NBASE
              XS(I,J)=PT(I,J)
  540     CONTINUE
          DO 550 J=1, NOBJ
              FXS(I,J)=IOBJ(I,J)
  550     CONTINUE
          IF (IPRINT) THEN
              WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3))")') NBASE,NOBJ
              WRITE(6,FMT) ((XS(I,J)),J=1,NBASE), (FXS(I,J),J=1,NOBJ)
          ENDIF
  560 CONTINUE
      RETURN
      END
C***********************************************************************
C
C     MULTI-OBJECTIVE GENETIC ALGORITHM
C
C***********************************************************************
      SUBROUTINE NSGAII(NGEN,POPSIZE,NBASE,NOBJ,XS,FXS)
C***********************************************************************
      DOUBLE PRECISION SEED0
      EXTERNAL FUN
      INTEGER MAXPOP, MAXBAS, MAXOBJ, NGEN, POPSIZE, NBASE, NOBJ
      PARAMETER (MAXPOP=2000, MAXBAS=200, MAXOBJ=3)
      INTEGER I, J, K, L, M, QSIZE
      INTEGER FSIZE(2*MAXPOP), NP(2*MAXPOP), RANK(2*MAXPOP)
      INTEGER Q(2*MAXPOP), SP(2*MAXPOP,2*MAXPOP)
      INTEGER SPSIZE(2*MAXPOP), F(2*MAXPOP,2*MAXPOP), ISORT(2*MAXPOP)
      INTEGER IMIN, IAUX, PAR1, PAR2, BETTER, WORSE
C      INTEGER COP
      REAL MAXXS, FX, UN, NO(1), CROSTYP, RATEMUT, TOTFIT, SUMFIT
C      REAL PCM
      REAL MAXFIT, XS(MAXPOP,MAXBAS), FXS(MAXPOP,MAXOBJ)
      REAL X(MAXBAS), PT(MAXPOP,MAXBAS), QT(MAXPOP,MAXBAS)
      REAL RT(2*MAXPOP,MAXBAS), IOBJ(2*MAXPOP,MAXOBJ), IDIST(2*MAXPOP)
      REAL ETA, BETA, QLIM, DL(MAXPOP), MAXID, SIGL, SUMM, SUM2, DELTAL
      CHARACTER (LEN=20) FILENAME
      INTEGER IDX, LDL, NDOM
      LOGICAL IPRINT, SNAP, CONVDEL, CONVPS, CONVDOM
      PARAMETER (ETA=8.0, LDL=40, DELTAL=0.02)
      CHARACTER (LEN=100) FMT
C      TYPE OF CROSSOVER (NONE:0; UNIFORM:1; FITNESS ORIENTED:2):
      PARAMETER (CROSTYP=1)
C      MUTATION RATE (WORKS BETTER BETWEEN 1%-2%):
      PARAMETER (RATEMUT=0.02)
C
C     TO PRINT RESULTS IPRINT=.TRUE., OTHERWISE, IPRINT=.FALSE.
C
      IPRINT=.FALSE.
C
C     F       SETS FI
C     FSIZE   SET FI SIZE
C     FX      THROUGHPUT
C     MAXXS   UPPER BOUND FOR XS VECTOR
C     IS      SET IS
C     NGEN    NUMBER OF GENERATIONS
C     POPSIZE POPULATION SIZE
C     MAXPOP  MAXIMUM POPULATION SIZE
C     NBASE   NUMBER OF BASES IN THE GENETIC CODE
C     NOBJ    NUMBER OF OBJECTIVE FUNCTIONS
C     MAXBAS  MAXIMUM NUMBER OF BASES
C     MAXOBJ  MAXIMUM NUMBER OF OBJECTIVES
C     PT      SET PT
C     QT      SET QT
C     Q       SET Q
C     QSIZE   SET Q SIZE
C     RANK    SET OF RANKS
C     NP      SET NP
C     SP      SET SP
C     SPSIZE  SET SP SIZE
C     RT      SET RT
C     IOBJ    OBJECTIVE FUNCTIONS OF SET RT
C     ISORT   SET SORTED I
C     IAUX    AUXILIARY VARIABLE FOR SORTING
C     IMIN    AUXILIARY VARIABLE FOR SORTING
C     SEED0   SEED FOR RANDOM NUMBER GENERATION
C     FUN     IMPLEMENTS THE GENERALIZED EXPANSION METHOD
C     XS      VECTOR OF CAPACITIES
C     FXS     VECTOR OF FUNCTION VALUES
C     UN      UNIFORMLY DISTRIBUTED RANDOM NUMBER
C     CROSTYP GA TYPE OF CROSSOVER (0, NO; 1, UNIFORM; 2, FITNESS ORIENTED)
C     RATEMUT GA MUTATION RATE (WORKS BETTER BETWEEN 1%-2%)
C
      IF (IPRINT) WRITE(6,*) 'NSGAII'
C
C     SET SEED FOR RANDOM NUMBER GENERATION
C
      SEED0=246813579.d0
C     SEED0=135792468.d0
C
C     GENERATE INITIAL POPULATION
C
      MAXXS=60
C      UN=1.0
      DO 20 I=1, POPSIZE
C          CALL UNI(SEED0,NO,1)
C          UN=NO(1)
          DO 10 J=1, NBASE/2
C
C             PT'S
C
              CALL UNI(SEED0,NO,1)
              UN=NO(1)
              PT(I,J)=CEILING(UN*MAXXS)
C             WRITE(6,*) 'CEILING(UN*MAXXS)', CEILING(UN*MAXXS)
              CALL UNI(SEED0,NO,1)
              UN=NO(1)
              PT(I,NBASE/2+J)=(UN*MAXXS)
C             WRITE(6,*) '(UN*MAXXS)', (UN*MAXXS)
C
C             QT'S
C
              CALL UNI(SEED0,NO,1)
              UN=NO(1)
              QT(I,J)=CEILING(UN*MAXXS)
C             WRITE(6,*) 'CEILING(UN*MAXXS)', CEILING(UN*MAXXS)
              CALL UNI(SEED0,NO,1)
              UN=NO(1)
              QT(I,NBASE/2+J)=(UN*MAXXS)
C             WRITE(6,*) '(UN*MAXXS)', (UN*MAXXS)
   10     CONTINUE
          IF (IPRINT) THEN
              WRITE(FMT,'("(A,I0,A,",I0,"(F10.4))")') NBASE
              WRITE(6,FMT) ' PT(',I,'): ',(PT(I,J),J=1,NBASE)
              WRITE(FMT,'("(A,I0,A,",I0,"(F10.4))")') NBASE
              WRITE(6,FMT) ' QT(',I,'): ',(QT(I,J),J=1,NBASE)
          ENDIF
   20 CONTINUE
C
C     GENERAL INITIALIZATIONS
C
      IDX=47
      CONVDEL=.TRUE.
      CONVPS=.TRUE.
      CONVDOM=.TRUE.
      DO 30 I=1, LDL
          DL(I)=0.0
   30 CONTINUE
C
C     ITERATE FOR NGEN GENERATIONS
C
      IF (IPRINT) WRITE(6,*) "      GEN      BEST_SOLUTION"
      I=0
   40 CONTINUE
C         IF (I.GE.2176.AND.I.LE.2180) THEN
C             IPRINT=1
C         ELSE
C             IPRINT=0
C         ENDIF
          IF (IPRINT) WRITE(6,*) 'NSGAII: GEN', I
C         WRITE(6,*) 'NSGAII: GEN', I
C
C         DEFINE SNAPSHOTS
C
          SNAP=(I.EQ.0.OR.I.EQ.10.OR.I.EQ.100.OR.I.EQ.NGEN)
C
C         COMBINE PARENT AND OFFSPRING POPULATION
C
          DO 70 J=1, POPSIZE
              DO 60 K=1, NBASE
C                 RT(2*J-1,K)=PT(J,K)
C                 RT(2*J,K)=QT(J,K)
                  RT(J,K)=PT(J,K)
                  RT(J+POPSIZE,K)=QT(J,K)
   60         CONTINUE
   70     CONTINUE
C
C         COMPUTE OBJECTIVE FUNCTIONS
C
          IF (IPRINT) WRITE(6,*) 'RT  IOBJ'
          DO 90 J=1,2*POPSIZE
              DO 80 K=1, NBASE/2
                  X(K)=RT(J,K)
                  X(NBASE/2+K)=RT(J,NBASE/2+K)
   80         CONTINUE
C             WRITE(6,*) 'NSGAII: CALLING FUN(RT,FX):'
              CALL FUN(X,FX,IPRINT)
C             UPDATE INPUT VECTOR
              DO 82 K=1, NBASE/2
                  RT(J,K)=X(K)
                  RT(J,NBASE/2+K)=X(NBASE/2+K)
   82         CONTINUE
C             COMPUTING OBJECTIVE FUNCTIONS
              IOBJ(J,1)=0.0
              IOBJ(J,2)=0.0
              DO 84 K=1, NBASE/2
                  IOBJ(J,1)=IOBJ(J,1)+X(K)
                  IOBJ(J,2)=IOBJ(J,2)+X(NBASE/2+K)
   84         CONTINUE
              IOBJ(J,3)=FX
              IF (IPRINT) THEN
               WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3))")') NBASE,NOBJ
               WRITE(6,FMT) (RT(J,K),K=1,NBASE),(IOBJ(J,K),K=1,NOBJ)
              ENDIF
   90     CONTINUE
C
C         PERFORM FAST NONDOMINATED SORTING
C
          DO 100 J=1, 2*POPSIZE
              FSIZE(J)=0
  100     CONTINUE
          DO 140 J=1, 2*POPSIZE
              SPSIZE(J)=0
              NP(J)=0
              DO 130 K=1, 2*POPSIZE
C                 DOES J DOMINATE K?
                  BETTER=0
                  WORSE=0
                  DO 110 L=1, NOBJ
C                     WRITE(6,*) 'IOBJ(J), IOBJ(K)', IOBJ(J,L),IOBJ(K,L)
                      IF (IOBJ(J,L).LT.IOBJ(K,L)) THEN
                          BETTER=BETTER+1
                      ELSE IF (IOBJ(J,L).GT.IOBJ(K,L)) THEN
                          WORSE=WORSE+1
                      ENDIF
  110             CONTINUE
C                 WRITE(6,*) 'BETTER', BETTER, ' WORSE', WORSE
                  IF (BETTER.GT.0.AND.WORSE.EQ.0) THEN
C                     WRITE(6,*) J,' DOMINATES', K, '->UPDATE SP'
                      SPSIZE(J)=SPSIZE(J)+1
                      SP(J,SPSIZE(J))=K
                  ELSE
C                     DOES K DOMINATE J?
                      BETTER=0
                      WORSE=0
                      DO 120 L=1, NOBJ
C                         WRITE(6,*) 'IOBJ(K), IOBJ(J)', IOBJ(K,L),IOBJ(J,L)
                          IF (IOBJ(K,L).LT.IOBJ(J,L)) THEN
                              BETTER=BETTER+1
                          ELSE IF (IOBJ(K,L).GT.IOBJ(J,L)) THEN
                              WORSE=WORSE+1
                          ENDIF
  120                 CONTINUE
C                     WRITE(6,*) 'BETTER', BETTER, ' WORSE', WORSE
                      IF ((BETTER.GT.0).AND.(WORSE.EQ.0)) THEN
C                         WRITE(6,*) K,' DOMINATES', J, '->UPDATE NP'
                          NP(J)=NP(J)+1
                      ENDIF
                  ENDIF
  130         CONTINUE
              IF (IPRINT) THEN
                  IF (SPSIZE(J).GT.0) THEN
                   WRITE(FMT,'("(A,I0,A,",I0,"I4)")') SPSIZE(J)
                   WRITE(6,FMT)' SP(',J, '): ',((SP(J,K)),K=1,SPSIZE(J))
C                   WRITE(6,132) J, ((SP(J,K)),K=1,SPSIZE(J))
C  132              FORMAT(' SP(',I0,'): ', <SPSIZE(J)>I4)
                  ENDIF
              ENDIF
              IF (NP(J).EQ.0) THEN
                  RANK(J)=1
                  FSIZE(1)=FSIZE(1)+1
                  F(1,FSIZE(1))=J
              ENDIF
  140     CONTINUE
C
          IF (IPRINT) THEN
              WRITE(6,*) 'NP: '
              WRITE(6,*) (NP(J),J=1,2*POPSIZE)
              WRITE(FMT,'("(A,",I0,"I4)")') FSIZE(1)
              WRITE(6,FMT) ' F(1): ', (F(1,K),K=1,FSIZE(1))
C              WRITE(6,142) (F(1,K),K=1,FSIZE(1))
C 142          FORMAT(' F(1): ', <FSIZE(1)>I4)
          ENDIF
C
          J=1
  150     CONTINUE
              QSIZE=0
C             WRITE(6,*)'P',(F(J,K),K=1,FSIZE(J))
              DO 160 K=1, FSIZE(J)
C                 P=F(J,K)
C                 WRITE(6,*)'Q',(SP(F(J,K),L),L=1,SPSIZE(F(J,K)))
                  DO 154 L=1, SPSIZE(F(J,K))
C                     Q=SP(P,L)
C                     Q=SP(F(J,K),L)
                      NP(SP(F(J,K),L))=NP(SP(F(J,K),L))-1
C
C                     ATTENTION: IMPORTANT CHANGE HERE, O TO 0
C
                      IF (NP(SP(F(J,K),L)).EQ.0) THEN
                          RANK(SP(F(J,K),L))=J+1
                          QSIZE=QSIZE+1
                          Q(QSIZE)=SP(F(J,K),L)
                      ENDIF
  154             CONTINUE
C                 WRITE(6,*) 'NP', (NP(L),L=1,2*POPSIZE)
  160         CONTINUE
              J=J+1
              FSIZE(J)=0
              DO 170 K=1, QSIZE
                  FSIZE(J)=FSIZE(J)+1
                  F(J,FSIZE(J))=Q(K)
  170         CONTINUE
              IF (IPRINT) THEN
                  IF (FSIZE(J).GT.0) THEN
                    WRITE(FMT,'("(A,I0,A,",I0,"I4)")') FSIZE(J)
                    WRITE(6,FMT) ' F(',J,'): ', (F(J,K),K=1,FSIZE(J))
C                    WRITE(6,172) J, ((F(J,K)),K=1,FSIZE(J))
C 172                FORMAT(' F(',I0,'): ', <FSIZE(J)>I4)
                  ENDIF
              ENDIF
          IF ((J.LT.(2*POPSIZE)).AND.(FSIZE(J).NE.0)) GOTO 150
C
C         COMPUTE CROWDING DISTANCES
C
          IF (IPRINT) WRITE(6,*) 'IDIST:'
          DO 180 J=1, 2*POPSIZE
              IDIST(J)=0.0
  180     CONTINUE
          DO 230 J=1, NOBJ
C             SORT BY J-TH OBJECTIVE
C             WRITE(6,*) 'OBJ', J
              DO 190 K=1, 2*POPSIZE
                  ISORT(K)=K
  190         CONTINUE
C             WRITE(6,*) 'IOBJ', (IOBJ(ISORT(L),J),L=1,2*POPSIZE)
              DO 210 K=1, 2*POPSIZE
                  IMIN=K
                  DO 200 L=K+1, 2*POPSIZE
                      IF (IOBJ(ISORT(L),J).LT.IOBJ(ISORT(IMIN),J)) THEN
                          IMIN=L
C                     ELSE
C                         IF (IOBJ(ISORT(L),J).EQ.IOBJ(ISORT(IMIN),J)
C     +                       .AND.IDIST(ISORT(L))
C     +                       .GT.IDIST(ISORT(IMIN))) THEN
C                             IMIN=L
C                         ENDIF
                      ENDIF
  200             CONTINUE
                  IAUX=ISORT(K)
                  ISORT(K)=ISORT(IMIN)
                  ISORT(IMIN)=IAUX
  210         CONTINUE
C
C             PRINT SNAPSHOT
C
C             IF (SNAP) THEN
C                 WRITE(6,*) 'NSGAII: SORTED IOBJ AFTER GENERATION', I
C                 WRITE(6,*) (IOBJ(ISORT(K),J),K=1,2*POPSIZE)
C             ENDIF
C             UPDATE CROWDING DISTANCES
              IDIST(ISORT(1))=1E+38
              IDIST(ISORT(2*POPSIZE))=1E+38
C             WRITE(6,*) 'FMAX', IOBJ(ISORT(1),J)
C             WRITE(6,*) 'FMIN', IOBJ(ISORT(2*POPSIZE),J)
              DO 220 K=2, (2*POPSIZE)-1
                  IDIST(ISORT(K))=IDIST(ISORT(K))+
     +                (IOBJ(ISORT(K+1),J)-IOBJ(ISORT(K-1),J))/
     +                (IOBJ(ISORT(2*POPSIZE),J)-IOBJ(ISORT(1),J))
  220         CONTINUE
C
C             PRINT SNAPSHOT
C
C             IF (SNAP) THEN
C                 WRITE(6,*) 'NSGAII: IDIST AFTER GENERATION', I
C                 WRITE(6,*) (IDIST(ISORT(K)), K=1,2*POPSIZE)
C             ENDIF
  230     CONTINUE
          IF (IPRINT) THEN
              WRITE(FMT,'("(",I0,"(E10.3))")') 2*POPSIZE
              WRITE(6,FMT) (IDIST(ISORT(K)), K=1,2*POPSIZE)
C              WRITE(6,232) (IDIST(ISORT(K)), K=1,2*POPSIZE)
C  232         FORMAT(<2*POPSIZE>(E10.3))
          ENDIF
C
C         PRINT SNAPSHOT OF POPULATION AT EACH LEVEL
C
          IF (SNAP) THEN
            IDX=IDX+1
            FILENAME='SNAP-'//CHAR(IDX)//'.TXT'
            OPEN(UNIT=7, FILE=FILENAME, STATUS='UNKNOWN')
C            WRITE(6,*) 'SNAPSHOT WRITTEN IN FILE ', FILENAME
            WRITE(7,*)
     +            'NSGAII: POPULATION AT EACH LEVEL AFTER GENERATION', I
            J=1
  240       CONTINUE
              DO 250 L=1, FSIZE(J)
                WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3),I4,E10.3)")')
     +                NBASE, NOBJ
                WRITE(7,FMT) ((RT(F(J,L),K)),K=1,NBASE),
     +            (IOBJ(F(J,L),K),K=1,NOBJ), J, IDIST(F(J,L))
  250         CONTINUE
              J=J+1
            IF (FSIZE(J).NE.0) GOTO 240
            CLOSE(UNIT=7)
          ENDIF
C
C         TRACK CONVERGENCE
C
C         COMPUTE MAXIMAL CROWDING DISTANCE
C
          MAXID=0.0
          DO 260 J=1, 2*POPSIZE
              IF (IDIST(J).GT.MAXID.AND.IDIST(J).LT.1E+38) THEN
                  MAXID=IDIST(J)
              ENDIF
  260     CONTINUE
C         WRITE(6,*) 'MAXID:', MAXID
C
C         STORE VALUE
C
          DL(MOD(I,LDL)+1)=MAXID
C         WRITE(6,*) 'POSITION:', MOD(I,LDL)+1
          IF (IPRINT) WRITE(6,*) 'DL:', (DL(J),J=1,LDL)
C
C         COMPUTE SIGL
C
          SUMM=0.0
          SUM2=0.0
          IF (I.GE.LDL) THEN
              DO 280 J=1, LDL
                  SUMM=SUMM+DL(J)
                  SUM2=SUM2+DL(J)**2
  280         CONTINUE
              SIGL=SQRT(ABS(SUM2-SUMM*SUMM/LDL)/LDL)
          ELSEIF (I.NE.0) THEN
              DO 270 J=1, I+1
                  SUMM=SUMM+DL(J)
                  SUM2=SUM2+DL(J)**2
  270         CONTINUE
              SIGL=SQRT(ABS(SUM2-SUMM*SUMM/(I+1))/(I+1))
          ELSE
              SIGL=0.0
          ENDIF
          IF (IPRINT) WRITE(6,*) 'SUM2, SUMM, SIGL:', SUM2, SUMM, SIGL
          IF (CONVDEL.AND.SIGL.LT.DELTAL.AND.I.GE.LDL) THEN
            CONVDEL=.FALSE.
C           WRITE(6,*)'NSGAII: DELTA_L BELOW THRESHOLD AT GENERATION', I
          ENDIF
C
C         USE SELECTION, CROSSOVER, AND MUTATION TO CREATE NEW POPULATION
C
C         SELECTION
C
          J=1
          K=0
          IF (IPRINT) WRITE(6,*) 'PT_INI:'
  290     IF (K+(FSIZE(J)).GT.(POPSIZE)) GOTO 320
              IF (IPRINT) WRITE(6,*) ' FSIZE(',J,')=', FSIZE(J)
              DO 310 L=1, FSIZE(J)
                  K=K+1
                  DO 300 M=1, NBASE
                      PT(K,M)=RT(F(J,L),M)
  300             CONTINUE
                  IF (IPRINT) THEN
                      WRITE(FMT,'("(",I0,"(F10.4))")') NBASE
                      WRITE(6,FMT) (PT(K,M),M=1,NBASE)
                  ENDIF
  310         CONTINUE
              J=J+1
              GOTO 290
  320     CONTINUE
          IF (IPRINT) WRITE(6,*) 'PT    '
  330     IF (K.GE.POPSIZE) GOTO 360
C             WRITE(6,*) ' FSIZE(',J,')=', FSIZE(J)
C
C             FIND LARGEST CROWDING DISTANCE
C
C             ERROR DETECTED BY MEANS OF TOM'S EXAMPLE
C             IMAX=F(J,1)
C
              IMAX=1
              DO 340 L=2, FSIZE(J)
                  IF (IDIST(F(J,IMAX)).LT.IDIST(F(J,L))) THEN
                      IMAX=L
                  ENDIF
  340         CONTINUE
C             WRITE(6,*) 'IMAX', F(J,IMAX)
C
C             DELETE THIS CROWDING DISTANCE
C
              IDIST(F(J,IMAX))=-1.0
C
C             INCLUDE ELEMENT WITH LARGEST CROWDING DISTANCE
C
              K=K+1
              DO 350 M=1, NBASE
                  PT(K,M)=RT(F(J,IMAX),M)
  350         CONTINUE
              IF (IPRINT) THEN
                  WRITE(FMT,'("(",I0,"(F10.4))")') NBASE
                  WRITE(6,FMT) (PT(K,M),M=1,NBASE)
              ENDIF
              GOTO 330
  360     CONTINUE
C
C         FILL IT UP IF NECESSARY
C
          IF (IPRINT) WRITE(6,*) 'PT_COM:'
  370     IF (K.GE.POPSIZE) GOTO 390
              K=K+1
              CALL UNI(SEED0,NO,1)
              UN=NO(1)
              PAR1=CEILING(UN*2*POPSIZE)
              DO 380 M=1, NBASE
                  PT(K,M)=RT(PAR1,M)
  380         CONTINUE
              IF (IPRINT) THEN
                  WRITE(FMT,'("(",I0,"(F10.4))")') NBASE
                  WRITE(6,FMT) (PT(K,M),M=1,NBASE)
              ENDIF
              GOTO 370
  390     CONTINUE
C
C         COMPUTE OBJECTIVE FUNCTIONS (REDUNDANT)
C
          IF (IPRINT) WRITE(6,*) 'PT IOBJ'
          MAXFIT=0.0
          DO 410 J=1,POPSIZE
              DO 400 K=1, NBASE/2
                  X(K)=PT(J,K)
                  X(NBASE/2+K)=PT(J,NBASE/2+K)
  400         CONTINUE
C             WRITE(6,*)
C             WRITE(6,*) 'NSGAII: CALLING FUN(PT,FX):'
              CALL FUN(X,FX,IPRINT)
C             UPDATE INPUT VECTOR
              DO 402 K=1, NBASE/2
                  PT(J,K)=X(K)
                  PT(J,NBASE/2+K)=X(NBASE/2+K)
  402         CONTINUE
C             COMPUTING OBJECTIVE FUNCTIONS
              IOBJ(J,1)=0.0
              IOBJ(J,2)=0.0
              DO 404 K=1, NBASE/2
                  IOBJ(J,1)=IOBJ(J,1)+X(K)
                  IOBJ(J,2)=IOBJ(J,2)+X(NBASE/2+K)
  404         CONTINUE
              IOBJ(J,3)=FX
              IF (-IOBJ(J,3).GT.MAXFIT) THEN
                  MAXFIT=-IOBJ(J,3)
              ENDIF
              IF (IPRINT) THEN
                  WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3))")')
     +                  NBASE,NOBJ
                  WRITE(6,FMT) (PT(J,K),K=1,NBASE),(IOBJ(J,K),K=1,NOBJ)
              ENDIF
C             WRITE(6,406) (PT(J,K),K=1,NBASE),(IOBJ(J,K),K=1,NOBJ)
 410      CONTINUE
C
C         TRACK CONVERGENCE
C
C         COMPUTE NUMBER OF DOMINATED SOLUTIONS (AT PARETO SET ONLY)
C
          NDOM=0
          DO 440 K=1, MIN(POPSIZE,FSIZE(1))
C             DOES FXS DOMINATES K?
              BETTER=0
              WORSE=0
              DO 420 L=1, NOBJ
C                 WRITE(6,*) 'FXS, IOBJ(K)', FXS(1,L),IOBJ(K,L)
                  IF (FXS(1,L).LT.IOBJ(K,L)) THEN
                      BETTER=BETTER+1
                  ELSE IF (FXS(1,L).GT.IOBJ(K,L)) THEN
                      WORSE=WORSE+1
                  ENDIF
  420         CONTINUE
C             WRITE(6,*) 'BETTER', BETTER, ' WORSE', WORSE
              IF (BETTER.GT.0.AND.WORSE.EQ.0) THEN
C                 WRITE(6,*) 'FXS DOMINATES', K, '->UPDATE NDOM'
                  NDOM=NDOM+1
              ELSE
C                 DOES K DOMINATE FXS?
                  BETTER=0
                  WORSE=0
                  DO 430 L=1, NOBJ
C                     WRITE(6,*) 'IOBJ(K), FXS', IOBJ(K,L),FXS(1,L)
                      IF (IOBJ(K,L).LT.FXS(1,L)) THEN
                          BETTER=BETTER+1
                      ELSE IF (IOBJ(K,L).GT.FXS(1,L)) THEN
                          WORSE=WORSE+1
                      ENDIF
  430             CONTINUE
C                 WRITE(6,*) 'BETTER', BETTER, ' WORSE', WORSE
                  IF (CONVDOM.AND.BETTER.GT.0.AND.WORSE.EQ.0) THEN
                      CONVDOM=.FALSE.
                      WRITE(6,*)
     +                     'NSGAII: WARNING: FXS HAS BEEN DOMINATED'
                  ENDIF
              ENDIF
  440     CONTINUE
C
C         PRINT OUT CONVERGENCE RESULTS
C
C         WRITE(6,*) I, SIGL, NDOM
C
C         WARN WHEN CONVERGENCE IS FIRST REACHED
C
          IF (CONVPS.AND.NDOM.EQ.0) THEN
            CONVPS=.FALSE.
C           WRITE(6,*)'NSGAII: TARGET IN PARETO SET AT GENERATION', I
          ENDIF
C
C         PRINT SNAPSHOT OF THE RESULT (PARETO-SET ONLY) AFTER GENERATION I
C
          IF (SNAP) THEN
              IDX=IDX+1
              FILENAME='SNAP-'//CHAR(IDX)//'.TXT'
              OPEN(UNIT=7, FILE=FILENAME, STATUS='UNKNOWN')
C              WRITE(6,*) 'SNAPSHOT WRITTEN IN FILE ', FILENAME
              WRITE(7,*) 'NSGAII: RESULT AFTER GENERATION', I
              DO 450 J=1, MIN(POPSIZE,FSIZE(1))
                WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3))")')NBASE,NOBJ
                WRITE(7,FMT) ((PT(J,K)),K=1,NBASE),(IOBJ(J,K),K=1,NOBJ)
  450         CONTINUE
              CLOSE(UNIT=7)
          ENDIF
C
C         CROSSOVER AND MUTATION
C
          IF (IPRINT) WRITE(6,*) 'QT'
          DO 510 J=1, POPSIZE
C             NO CROSSOVER
              IF (CROSTYP.EQ.0) THEN
C                 WRITE(6,*) 'NO CROSSOVER'
                  PAR1=J
                  PAR2=J
C             UNIFORMLY DISTRIBUTED CROSSOVER
              ELSEIF (CROSTYP.EQ.1) THEN
C                 WRITE(6,*) 'UNIFORMLY DISTRIBUTED CROSSOVER'
                  CALL UNI(SEED0,NO,1)
                  UN=NO(1)
                  PAR1=CEILING(UN*POPSIZE)
  460             CONTINUE
                      CALL UNI(SEED0,NO,1)
                      UN=NO(1)
                      PAR2=CEILING(UN*POPSIZE)
                  IF (PAR1.EQ.PAR2) GOTO 460
C             FITNESS ORIENTED CROSSOVER (ONLY FOR SINGLE-OBJECTIVE)
              ELSEIF (CROSTYP.EQ.2) THEN
C                 WRITE(6,*) 'FITNESS ORIENTED CROSSOVER'
                  TOTFIT=0.0
                  DO 470 K=1,POPSIZE
                      TOTFIT=TOTFIT+(1.1*MAXFIT-(-IOBJ(K,2)))
  470             CONTINUE
C                 CHOOSE PARENT 1
                  CALL UNI(SEED0,NO,1)
                  UN=NO(1)
                  PAR1=0
                  SUMFIT=0.0
C                 WRITE(6,*) 'TOTFIT*UN', TOTFIT*UN
  480             CONTINUE
                      PAR1=PAR1+1
                      SUMFIT=SUMFIT+(1.1*MAXFIT-(-IOBJ(PAR1,2)))
C                     WRITE(6,*) 'PAR1', PAR1, ' SUMFIT', SUMFIT
                  IF (PAR1+1.LT.POPSIZE.AND.
     +                SUMFIT+(1.1*MAXFIT-(-IOBJ(PAR1+1,2))).LE.
     +                TOTFIT*UN) GOTO 480
                  IF (SUMFIT.LT.TOTFIT*UN) THEN
                      PAR1=PAR1+1
                      SUMFIT=SUMFIT+(1.1*MAXFIT-(-IOBJ(PAR1,2)))
C                     WRITE(6,*) 'PAR1', PAR1, ' SUMFIT', SUMFIT
                  ENDIF
C                 CHOOSE PARENT 2
                  CALL UNI(SEED0,NO,1)
                  UN=NO(1)
                  PAR2=0
                  SUMFIT=0.0
C                 WRITE(6,*) 'TOTFIT*UN', TOTFIT*UN
  490             CONTINUE
                      PAR2=PAR2+1
                      SUMFIT=SUMFIT+(1.1*MAXFIT-(-IOBJ(PAR2,2)))
C                     WRITE(6,*) 'PAR2', PAR2, ' SUMFIT', SUMFIT
                  IF (PAR2+1.LT.POPSIZE.AND.
     +                SUMFIT+(1.1*MAXFIT-(-IOBJ(PAR2+1,2))).LE.
     +                TOTFIT*UN) GOTO 490
                  IF (SUMFIT.LE.TOTFIT*UN) THEN
                      PAR2=PAR2+1
                      SUMFIT=SUMFIT+(-IOBJ(PAR2,2))
C                     WRITE(6,*) 'PAR2', PAR2, ' SUMFIT', SUMFIT
                  ENDIF
              ELSE
                  WRITE(6,*) 'NSGAII: ERROR: UNKNOWN CROSSOVER'
                  STOP
              ENDIF
C              WRITE(6,*) 'PAR1    ', PAR1, '  PAR2    ', PAR2
C              WRITE(6,492) (PT(PAR1,K),K=1,NBASE)
C              WRITE(6,494) (PT(PAR2,K),K=1,NBASE)
C  492         FORMAT(' PT(PAR1):', <NBASE>(F,X))
C  494         FORMAT(' PT(PAR2):', <NBASE>(F,X))
C
C             SELECT CROSSOVER POSITION (LUIZ DUCZMAL)
C
C             CALL UNI(SEED0,UN,1)
C             COP=CEILING(UN*NBASE)
C             WRITE(6,*) 'COP', COP
C
C             CROSSOVER & MUTATION
C
              DO 500 K=1, NBASE
C
C                 CROSSOVER TYPE #1 (LUIZ DUCZMAL)
C
C                 PCM=MOD((K-COP+NBASE),NBASE)
C                 QT(J,K)=(PCM*PT(PAR1,K)+(NBASE-1-PCM)*PT(PAR2,K))/
C     +               FLOAT(NBASE-1)
C                 WRITE(6,*) 'PCM', K, ':', PCM
C                 WRITE(6,*) QT(J,K), '=',
C     +               PCM/FLOAT(NBASE-1), '*', PT(PAR1,K), '+',
C     +               (NBASE-1-PCM)/FLOAT(NBASE-1), '*', PT(PAR2,K)
C
C                 CROSSOVER TYPE #2 (GRAHAN KENDALL)
C
C                 QT(J,K)=0.5*PT(PAR1,K) + 0.5*PT(PAR2,K)
C                 WRITE(6,*) QT(J,K), '= .5*', PT(PAR1,K), '+ .5*', PT(PAR2,K)
C
C                 SIMULATED BINARY CROSSOVER, 0.5*100% (DEB & AGRAWAL, 1995)
C
                  CALL UNI(SEED0,NO,1)
                  UN=NO(1)
                  IF (UN.LT.0.5) THEN
                      CALL UNI(SEED0,NO,1)
                      UN=NO(1)
                      IF (UN.LE.0.5) THEN
                          BETA=(2*UN)**(1/(ETA+1))
                      ELSEIF (UN.LT.(1-1E-08)) THEN
                          BETA=(1/(2*(1-UN)))**(1/(ETA+1))
                      ELSE
                          BETA=(1/(2*(1E-08)))**(1/(ETA+1))
                      ENDIF
                      IF (IPRINT) WRITE(6,*) 'UN:', UN, 'BETA:', BETA
                      QT(J,K)=0.5*((1+BETA)*PT(PAR1,K)+
     +                            (1-BETA)*PT(PAR2,K))
                  ELSE
                      QT(J,K)=PT(PAR1,K)
                  ENDIF
C
C                 MUTATE, RATEMUT*100%
C
                  CALL UNI(SEED0,NO,1)
                  UN=NO(1)
                  IF (UN.LT.RATEMUT) THEN
C
C                     DUCZMAL IDEA
C
C                     CALL UNI(SEED0,UN,1)
C                     QT(J,K)=QT(J,K)+(-1)**CEILING(2.*UN)
C                     WRITE(6,*) 'MUTATION:'
C                     WRITE(6,*) '(-1)**CEILING(2.*UN)', (-1)**CEILING(2.*UN)
C
C                     NORMAL MUTATION (DEB & AGRAWAL, 1995)
C
                      CALL NORM(SEED0,NO,1)
                      UN=NO(1)
                      QT(J,K)=QT(J,K)+UN
C                     WRITE(6,*) 'MUTATION: K, UN:', K, UN
                  ENDIF
C
C                 ADJUSTMENT TO INTEGRALITY, NON-NEGATIVITY, AND OVERFLOW
C
                  IF (K.LE.NBASE/2) THEN
                      QT(J,K)=CEILING(QT(J,K))
                      QLIM=100.0
                      IF (QT(J,K).GT.QLIM) THEN
                          QT(J,K)=QLIM-ABS(QT(J,K)-QLIM)
C                         WRITE(6,*) K, '-TH CORRECTED TO:', QT(J,K)
                      ENDIF
                      QLIM=1.0
                      IF (QT(J,K).LT.QLIM) THEN
                          QT(J,K)=QLIM+ABS(QT(J,K)-QLIM)
C                         WRITE(6,*) K, '-TH CORRECTED TO:', QT(J,K)
                      ENDIF
                  ELSE
                      QLIM=50.0
                      IF (QT(J,K).GT.QLIM) THEN
                          QT(J,K)=QLIM-ABS(QT(J,K)-QLIM)
C                         WRITE(6,*) K, '-TH CORRECTED TO:', QT(J,K)
                      ENDIF
                      QLIM=0.0
                      IF (QT(J,K).LT.QLIM) THEN
                          QT(J,K)=QLIM+ABS(QT(J,K)-QLIM)
C                         WRITE(6,*) K, '-TH CORRECTED TO:', QT(J,K)
                      ENDIF
                  ENDIF
  500         CONTINUE
C
C             RESULT
C
              IF (IPRINT) THEN
                  WRITE(FMT,'("(A,I0,A,",I0,"(F10.4))")') NBASE
                  WRITE(6,FMT) ' QT(',J,'):  ',((QT(J,K)),K=1,NBASE)
              ENDIF
  510     CONTINUE
          I=I+1
      IF (I.LE.NGEN) GOTO 40
C  520 CONTINUE
C
C     FINAL RESULT
C
      IF (IPRINT) THEN
          WRITE(6,*) 'NSGAII RESULT AFTER'
          WRITE(6,*) 'GENERATIONS', NGEN
          WRITE(6,*) 'POPULATION', POPSIZE
          WRITE(6,*) 'XS FXS'
      ENDIF
      DO 550 I=1, POPSIZE
          DO 530 J=1, NBASE
              XS(I,J)=PT(I,J)
  530     CONTINUE
          DO 540 J=1, NOBJ
              FXS(I,J)=IOBJ(I,J)
  540     CONTINUE
          IF (IPRINT) THEN
              WRITE(FMT,'("(",I0,"(F10.4),",I0,"(E10.3))")') NBASE,NOBJ
              WRITE(6,FMT) ((XS(I,J)),J=1,NBASE), (FXS(I,J),J=1,NOBJ)
          ENDIF
  550 CONTINUE
      RETURN
      END
C***********************************************************************
C
C      UNI GENERATES A VECTOR R OF N UNIFORM(0,1) OBSERVATIONS
C      VARIATION ON IMSL ROUTINE GGUBS, COURTESY MARTIN TANNER
C
C***********************************************************************
      SUBROUTINE UNI(DSEED,R,N)
C***********************************************************************
      REAL R(1)
      DOUBLE PRECISION DSEED,D2P31M,D2P31
      DATA D2P31M,D2P31 /2147483647.d0,2147483648.d0/
C
      DO 10 I=1,N
            DSEED=DMOD(16807.D0*DSEED,D2P31M)
            R(I)=REAL(DSEED/D2P31)
   10 CONTINUE
C
      RETURN
      END
C***********************************************************************
C
C      NORM GENERATES A VECTOR R OF N NORMAL(0,1) OBSERVATIONS
C
C***********************************************************************
      SUBROUTINE NORM(DSEED,R,N)
C***********************************************************************
      REAL R(1), U(2), Y(2)
      DOUBLE PRECISION DSEED
C
      I=0
   10 CONTINUE
            CALL UNI(DSEED,U,2)
            Y(1)=-LOG(U(1))
            Y(2)=-LOG(U(2))
            IF ((Y(2)-(Y(1)-1)*(Y(1)-1)/2).GT.0) THEN
                  I=I+1
                  CALL UNI(DSEED,U,1)
                  IF (U(1).GT.0.5) THEN
                        R(I)=-Y(1)
                  ELSE
                        R(I)=Y(1)
                  ENDIF
            ENDIF
      IF (I.LT.N) GOTO 10
C
      RETURN
      END
C***********************************************************************
C
C      FUNCTION TO CALCULATE PERFORMANCE MEASURES (GEM)
C
C***********************************************************************
      SUBROUTINE FUN(W,FX,IPRT)
C***********************************************************************
      COMMON/ONE/XLEFF,RATE(40),RH,RIC(40),RHO,XL(40),NQ,AS(40),ARH(40),
     + NARCS,NS(780),NF(780),STEST(40),SLAM(40),RP(780),ETEST(40),NRUN,
     + BETA(40),GAM(40),M,PC(40),PZ(40),EXL,TCAP,AMU,THTPUT,TOT,PENALTY
      COMMON/TWO/X(40),Y(40),
     + S(40),FY,N,KOUNTS,LIN,NDRV,DIRECT(40,40),BEFORE(40),
     + FIRST(40),SECND(40)
      COMMON/I/I
      COMMON/FOUR/WW(0:50)
      INTEGER NE,ITMAX,NTEST(50),COUNT,NQ,J,XBUF
      REAL R(0:40),FLOW(40,40),XGUESS(12),W(40),FX,TCAP,ROW,RMIN,ERRREL
      LOGICAL IPRT,IPRINT
      EXTERNAL FCN
      CHARACTER (LEN=100) FMT
C
C      TO PRINT RESULTS IPRINT=.TRUE., OTHERWISE, IPRINT=.FALSE.
C
      IPRINT=IPRT
      IPRINT=.FALSE.
C
C      EXL     EXPECTED QUEUE LENGTH
C      ITMAX   NUMBER OF ITERATIONS IN BROWN'S METHOD
C      M       FINITE QUEUE BUFFER SIZE
C      NARCS   NUMBER OF ARCS
C      NE      NUMBER OF SIMULTANEOUS EQUATIONS
C      NQ      NUMBER OF FINITE QUEUES
C      NS,NF   START AND END NODES FOR THE ARCS
C      NSNOD   NUMBER OF STARTING NODES
C      PC      PROBABILITY QUEUE IS AT CAPACITY
C      PZ      PROBABILITY OF 0 BUSY SERVERS AT A QUEUE
C      R       RHO AT NODE I
C      RATE    SERVICE RATE AT NODE I
C      RH      OVERALL SERVICE RATE AT THE HOLDING NODE
C      RIC     CAPACITY OF QUEUE I
C      RP      ROUTING PROBABILITIES ALONG ARCS
C      SLAM    INITIAL ARRIVAL RATE TO NODE I
C      SOUR    ARRAY THAT HOLDS THE SOURCE NODES
C      XGUESS  INITIAL GUESS FOR BROWN'S METHOD
C      XLEFF   EFFECTIVE ARRIVAL RATE (LAMBDA EFF)
C      FMT     FORMAT SPECIFIER
C
      IF (IPRINT) THEN
          WRITE(FMT,'("(A,",I0,"(F10.4))")') 2*NQ
          WRITE(6,FMT) ' FUN: W: ', (W(I),I=1,2*NQ)
C          WRITE(6,2) (W(I),I=1,2*NQ)
C    2     FORMAT(' FUN: W:', <2*NQ>(F,X))
          WRITE(6,*) 'FUN: NODE     XLEFF'
      ENDIF
C
C      INITIALIZATION OF PARAMETERS
C
      NE = 8
      ITMAX=100
      ERRREL=1.E-3
      TOT=0.0
      DO 10 I=1,NQ
            IF (W(I).LT.1) THEN
C                  PROBABLY A MISTAKE; ARBITRATE FX (HIGHEST VALUE) AND GET OUT
                  FX=NQ
                  PRINT*,'FX GOES TO', FX
                  RETURN
            ENDIF
            WW(I)=W(I)
            RATE(I)=W(NQ+I)
   10 CONTINUE
C
      XGUESS(1)=0.50
      XGUESS(2)=1.00
      XGUESS(3)=500.
      XGUESS(4)=1.00
      XGUESS(5)=0.01
      XGUESS(6)=2.00
      XGUESS(7)=0.01
      XGUESS(8)=0.05

C
C      LABEL ALL NODES AS UNEVALUATED
C
      DO 20 I=1,NQ
            NTEST(I)=0
   20 CONTINUE
C
C      IF ALL PREVIOUS NODES HAVING ARCS CONNECTED TO PRESENT NODE
C      HAVE BEEN EVALUATED, PRESENT NODE CAN BE EVALUATED AND FLOWS
C      ALONG INPUT ARCS ARE ADDED, AND BUFFER AUGMENTED.
C      OTHERWISE, GO TO NEXT NODE.
C      ALSO, IF A SOURCE NODE, FLOW IN EQUALS SLAM(I).
C
   30 DO 110 I=1,NQ
            IF (NTEST(I).EQ.1) GO TO 110
            IF (STEST(I).EQ.1) GO TO 50
            XLAM=0
            XBUF=0
            DO 40 J=1,NARCS
                  IF (NF(J).EQ.I.AND.NTEST(NS(J)).EQ.0) GO TO 110
                  IF (NF(J).EQ.I) THEN
C                        WRITE(6,*) 'FLOW', FLOW(NS(J),I)
                        XLAM=XLAM+FLOW(NS(J),I)
C                        WRITE(6,*) 'RIC(NS(J))', RIC(NS(J))
                        XBUF=XBUF+INT(RIC(NS(J)))
                  ENDIF
   40       CONTINUE
            XLEFF=XLAM
   50       IF (STEST(I).EQ.1) THEN
                  XLEFF=SLAM(I)
            ENDIF
C
C            SET BUFFER SIZE (INCLUDING THOSE IN SERVICE)
C
            RIC(I)=WW(I)
C
C            DETERMINE SERVICE RATE AT HOLDING NODE
C            MODIFIED ON JANUARY 22, 2002 FOR SPLIT AND MERGE TOPOLOGIES
C
C            HOLDING RATE ARH WAS PASSED TO PROGRAM 1/14/02
C            TO REFLECT HOLDING RATE AS A FUNCTION OF SQ. COEFF. OF VARIATION
C
            RH=ARH(I)
C
C            SET TOTAL CAPACITY SIZE (INCLUDING THOSE IN SERVICE)
C
            M=INT(WW(I))+1
C            WRITE(6,*) 'WW(I) ', WW(I)
C
C            FINITE BUFFER SIZE IS AUGMENTED FOR BLOCKING
C
C            ADJUST=WW(I-1)/2
C            M=WW(I)+MAX(CEILING(ADJUST),1)
C            PRINT *, 'M', M
C
C            FIND QUEUE LENGTH, BLOCKING PROBABILITY, THROUGHPUT
C
C            IF NECESSARY ADJUST ROW
C
            ROW=XLEFF/RATE(I)
            IF (ROW.GE.1.) THEN
C                  WRITE(6,*) 'FUN: RATE(I) IS UNFEASIBLE:', RATE(I)
                  RMIN=XLEFF+ERRREL
                  RATE(I)=RMIN+ABS(RATE(I)-RMIN)
                  ROW=XLEFF/RATE(I)
C                  WRITE(6,*) 'FUN: RATE(I) ADJUSTED TO:  ', RATE(I)
            ENDIF
C
C            PROGRAM CALL TO IMSL ROUTINE
C
C            WRITE(6,*) 'FUN: CALLING NEQNF:'
C            CALL NEQNF (FCN, ERRREL, NE, ITMAX, XGUESS, X, FNORM)
C            WRITE(6,*) 'FUN: RESULT NEQNF:'
C            WRITE(FMT,'("(A,",I,"(F10.4))")') NE
C            WRITE(6,FMT) ' FUN: X:',(X(J),J=1,NE)
C
C            PROGRAM CALL TO ALTERNATIVE ROUTINE
C            WRITTEN ON FEBRUARY 19, 2020
C
C            WRITE(6,*) 'FUN: CALLING ALTNEQNF:'
            CALL ALTNEQNF (FCN, ERRREL, NE, ITMAX, XGUESS, X, FNORM)
C            WRITE(6,*) 'FUN: RESULT ALTNEQNF:'
C            WRITE(FMT,'("(A,",I,"(F10.4))")') NE
C            WRITE(6,FMT) ' FUN: X:',(X(J),J=1,NE)
C
C            HERE WE AUGMENT THE EFFECTIVE ARRIVAL RATE
C
C   60       CONTINUE
            XLAM=XLEFF
            RHO=XLEFF/RH
            R(0)=0
            R(I)=RHO
            XLEFF=(1.0-X(8))*XLEFF
C
            TPUT = XLEFF
            IF (I.EQ.1) X(5)=0.0
C            AUGMENT=X(5)*(1.0-X(1))
C            FACTOR=(1.0-X(8))**(R(I-1)+R(I))
C            IF (RHO.GE.1.0) AUGMENT = 0.0
C            TPUT = TPUT + AUGMENT*FACTOR
            XLEFF = TPUT
C
C            DETERMINE FLOWS ALONG OUTPUT ARCS FROM EACH NODE
C
            DO 70 J=1,NARCS
                  IF (NS(J).EQ.I) THEN
                        FLOW(I,NF(J))=TPUT*RP(J)
C                        PRINT*,'FLOWOUT',I,NF(J),FLOW(I,NF(J))
                  ENDIF
   70       CONTINUE
C
C            OUTPUT RESULTS AT EACH NODE TO TERMINAL
C
            IF (IPRINT) THEN
              WRITE(6,80) I, XLEFF
   80         FORMAT(6X,I4,F10.4)
            ENDIF
C            WRITE(6,*) 'EXPANDED NETWORK PERFORMANCE MEASURES'
C            WRITE(6,90)
C   90       FORMAT(2X, 'LAMDA', 4X, 'RHO', 1X, 'CAPACITY', 5X, 'P1',
C     +            4X,'QUEUE', 1X, 'THRUPUT')
C            WRITE(6,100) XLAM, RHO, RIC(I), X(8), EXL, TPUT
C  100       FORMAT(F7.4, F7.4, F9.1, F7.4, F9.4, F8.4)
C
C            THE NODE HAS BEEN EVALUATED
C
            PC(I) = X(8)
            COUNT=COUNT+1
            NTEST(I)=1
C
C            IF AN END NODE, ADD TO TOTAL THRUPUT
C
            IF (ETEST(I).EQ.1) THEN
                  TOT=TOT+TPUT
            ENDIF
C
  110 CONTINUE
C
C      IF ALL NODES HAVE BEEN EVALUATED,STOP.
C      OTHERWISE, GO THROUGH PROCESS AGAIN TO EVALUATE
C      UNEVALUATED NODES.
C
      IF (COUNT.LT.NQ) GOTO 30
C
C      OBJECTIVE FUNCTION GENERATION FOR SINGLE OBJECTIVE APPROACH
C
C      FX=0.0
C      DO 120 I = 1,NQ
C            FX=FX+W(I)
C  120 CONTINUE
C      PENALTY=TOT-THTPUT
C      WRITE(6,*) 'FUN: PENALTY', PENALTY
C      FX=FX-AMU*(PENALTY)
C
C      YET ANOTHER OBJECTIVE FUNCTION
C
C      FX=-REV*TOT+COST*TCAPP
C
C      OBJECTIVE FUNCTION GENERATION FOR MULTIOBJETIVE APPROACH: THETA (ERRREL/10)
C
C      FX=INT(TOT/(ERRREL/10)+0.5)*(ERRREL/10)
C
C      WRITE(6,*) FX, TOT, TOT/ERRREL+0.5,
C     +      INT(TOT/ERRREL+0.5), INT(TOT/ERRREL+0.5)*ERRREL
C
C      OBJECTIVE FUNCTION GENERATION FOR MULTIOBJETIVE APPROACH: SUM PK (ERRREL/10)
C
      FX=0.0
      DO 120 I=1,NQ
            FX=FX+PC(I)
  120 CONTINUE
      FX=INT(FX/(ERRREL/10)+0.5)*(ERRREL/10)
C
C      WRITE(6,*) 'FUN: SUM PK ', FX
C
C      UPDATE W
C
      DO 130 I=1,NQ
            W(NQ+I)=RATE(I)
  130 CONTINUE
C
      IF (IPRINT) THEN
          WRITE(6,*) 'FUN: TOT ', TOT
          WRITE(6,*) 'FUN: TPUT', TPUT
          WRITE(6,*) 'FUN: FX  ', FX
      ENDIF
C
      RETURN
      END
C***********************************************************************
C
C      FUNCTION TO BE SOLVED BY BROWN'S METHOD
C
C***********************************************************************
      SUBROUTINE FCN(X,F,NE)
C***********************************************************************
      COMMON/ONE/XLEFF,RATE(40),RH,RIC(40),RHO,XL(40),NQ,AS(40),ARH(40),
     + NARCS,NS(780),NF(780),STEST(40),SLAM(40),RP(780),ETEST(40),NRUN,
     + BETA(40),GAM(40),M,PC(40),PZ(40),EXL,TCAP,AMU,THTPUT,TOT,PENALTY
      COMMON/FOUR/WW(0:50)
      COMMON/I/I
      INTEGER NE, IM
      REAL X(NE),F(NE)
      REAL ROW,TOP1,BOT1,T1,TOP2,BOT2,T2
C
C      PRINT*,'FCN: X: ', (X(J),J=1,NE)
C
C      PARAMETERS FOR M/G/1/K FORMULA
C
C      PRINT*, 'XLEFF=',XLEFF,' RATE =',RATE(I)
      ROW=XLEFF/RATE(I)
C      PRINT*,'ROW ', ROW
      TOP1 = 2.0 + SQRT(ROW)*AS(I) - SQRT(ROW) + 2.0*(M-1)
C      PRINT*,'TOP1 ', TOP1
      BOT1 = 2.0 + SQRT(ROW)*AS(I) - SQRT(ROW)
C      PRINT*,'BOT1 ', BOT1
      T1 = TOP1/BOT1
C      PRINT*,'T1 ', T1
      TOP2 = 2.0 + SQRT(ROW)*AS(I) - SQRT(ROW) + (M-1)
C      PRINT*,'TOP2 ', TOP2
      BOT2 = BOT1
C      PRINT*,'BOT2 ', BOT2
      T2 = 2.0*TOP2/BOT2
C      PRINT*,'T2 ', T2
      PC(I)=(ROW**T1 * (-1.0 + ROW))/(ROW**T2 -1.0)
C      PRINT*,'PC(I) ', PC(I)
C
C      MAX M TO AVOID OVERFLOW IN PK'
C
      IF (X(6).GT.1.) THEN
            IM=INT(38*LOG(10.)/LOG(X(6)))
            IF (IM.GT.M) THEN
                  IM=M
C                  WRITE(6,*) 'FCN: M UNCHANGED:', IM
            ELSE
C                  WRITE(6,*) 'FCN: M REDUCED TO:', IM
            ENDIF
      ELSE
            IM=M
C            WRITE(6,*) 'FCN: M UNCHANGED:', IM
      ENDIF
C
C      SIMULTANEOUS EQUATIONS
C
      XN1=((X(6)**IM)-(X(7)**IM))-((X(6)**(IM-1))-(X(7)**(IM-1)))
C      PRINT*,'XN1 ', XN1
      XD1=((X(6)**(IM+1))-(X(7)**(IM+1)))-((X(6)**IM)-(X(7)**IM))
C      PRINT*,'XD1 ', XD1
      F(1)=X(1)*(2.-(X(2)*XN1)/(RH*XD1))-1.
C      PRINT*,'F(1) ', F(1)
      F(2)=X(2)-X(4)+X(5)*(1.-X(1))
C      PRINT*,'F(2) ', F(2)
      F(3)=(X(2)+2*RH)**2-(4*(RH*X(2)))-X(3)
C      PRINT*,'F(3) ', F(3)
      F(4)=(1.-X(8))*XLEFF-X(4)
C      PRINT*,'F(4) ', F(4)
      F(5)=XLEFF-X(4)-X(5)
C      PRINT*,'F(5) ', F(5)
      F(6)=(((X(2)+2.*RH)+(X(3)**0.5))/(2.*RH))-X(6)
C      PRINT*,'F(6) ', F(6)
      F(7)=(((X(2)+2.*RH)-(X(3)**0.5))/(2.*RH))-X(7)
C      PRINT*,'F(7) ', F(7)
      F(8)=PC(I)-X(8)
C      PRINT*,'F(8) ', F(8)
C
      RETURN
      END
C***********************************************************************
C
C      ALTERNATIVE IMPLEMENTATION TO NEQNF
C      WRITTEN ON FEBRUARY 19, 2020
C
C***********************************************************************
      SUBROUTINE ALTNEQNF (FCN, ERRREL, NE, ITMAX, XGUESS, X, FNORM)
C***********************************************************************
      COMMON/ONE/XLEFF,RATE(40),RH,RIC(40),RHO,XL(40),NQ,AS(40),ARH(40),
     + NARCS,NS(780),NF(780),STEST(40),SLAM(40),RP(780),ETEST(40),NRUN,
     + BETA(40),GAM(40),M,PC(40),PZ(40),EXL,TCAP,AMU,THTPUT,TOT,PENALTY
      COMMON/FOUR/WW(0:50)
      COMMON/I/I
      INTEGER NE,ITMAX,IM,IT,J
      REAL ERRREL,XGUESS(NE),X(NE),FNORM,F(NE),XF(NE)
      REAL ROW,TOP1,BOT1,T1,TOP2,BOT2,T2
      EXTERNAL FCN
      LOGICAL GOON
C      CHARACTER (LEN=100) FMT
C
C      INITILIZE
C
      DO 10 J = 1,NE
            X(J)=XGUESS(J)
   10 CONTINUE
C      WRITE(FMT,'("(A,",I0,"(F10.4))")') NE
C      WRITE(6,FMT) ' ALTNEQNF: X:', (X(J),J=1,NE)
C
C      PERFORM ITERATIONS
C
      IT=0
   20 CONTINUE
            IT=IT+1
C            WRITE(6,*) 'ALTNEQNF: ITER ',IT
C
C            PARAMETERS FOR M/G/1/K FORMULA
C
            ROW=XLEFF/RATE(I)
            TOP1 = 2.0 + SQRT(ROW)*AS(I) - SQRT(ROW) + 2.0*(M-1)
            BOT1 = 2.0 + SQRT(ROW)*AS(I) - SQRT(ROW)
            T1 = TOP1/BOT1
            TOP2 = 2.0 + SQRT(ROW)*AS(I) - SQRT(ROW) + (M-1)
            BOT2 = BOT1
            T2 = 2.0*TOP2/BOT2
            PC(I)=(ROW**T1 * (-1.0 + ROW))/(ROW**T2 -1.0)
C
C            MAX M TO AVOID OVERFLOW IN PK'
C
            IF (X(6).GT.1.) THEN
                  IM=INT(38*LOG(10.)/LOG(X(6)))
                  IF (IM.GT.M) THEN
                        IM=M
C                        WRITE(6,*) 'ALTNEQNF: M UNCHANGED:', IM
                  ELSE
C                        WRITE(6,*) 'ALTNEQNF: M REDUCED TO:', IM
                  ENDIF
            ELSE
                  IM=M
C                  WRITE(6,*) 'ALTNEQNF: M UNCHANGED:', IM
            ENDIF
C
C            RECURSIVE EQUATIONS
C
            XN1=((X(6)**IM)-(X(7)**IM))-((X(6)**(IM-1))-(X(7)**(IM-1)))
            XD1=((X(6)**(IM+1))-(X(7)**(IM+1)))-((X(6)**IM)-(X(7)**IM))
C            PK'
            XF(1)=1/(2.-(X(2)*XN1)/(RH*XD1))
C            LAMBDA
            XF(2)=X(4)-X(5)*(1.-X(1))
C            Z
            XF(3)=(X(2)+2*RH)**2-(4*(RH*X(2)))
C            LAMBDA_J
            XF(4)=(1.-X(8))*XLEFF
C            LAMBDA_H
            XF(5)=XLEFF-X(4)
C            R2
            XF(6)=(((X(2)+2.*RH)+(X(3)**0.5))/(2.*RH))
C            R1
            XF(7)=(((X(2)+2.*RH)-(X(3)**0.5))/(2.*RH))
C            PK
            XF(8)=PC(I)
C
C            WRITE(FMT,'("(A,",I0,"(F10.4))")') NE
C            WRITE(6,FMT) ' ALTNEQNF: XF: ', (XF(J),J=1,NE)
C
C            VERIFY CONVERGENCE
C
            GOON=.FALSE.
            DO 30 J = 1, NE
                  IF (ABS(XF(J)-X(J)).GE.ERRREL) GOON=.TRUE.
                  X(J)=XF(J)
   30       CONTINUE
C
C            STOP CRITERIA
C
      IF (GOON.AND.IT.LT.ITMAX) GOTO 20
C
C      COMPUTE FNORM
C
      CALL FCN(X,F,NE)
      FNORM=0.
      DO 40 J = 1,NE
            FNORM=FNORM+F(J)*F(J)
   40 CONTINUE
C      WRITE(6,*) 'ALTNEQNF: FNORM ', FNORM
C
      RETURN
      END
C***********************************************************************
