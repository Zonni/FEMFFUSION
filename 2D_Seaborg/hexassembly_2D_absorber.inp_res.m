
% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'May 19 2021 21:58:01' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 27])  = 'hexassembly_2D_absorber.inp' ;
WORKING_DIRECTORY         (idx, [1: 26])  = '/nfs/home/jgj/projects/UPV' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode01' ;
CPU_TYPE                  (idx, [1: 46])  = 'AMD Ryzen Threadripper 3970X 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 137367609.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct  4 21:34:37 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct  4 21:35:17 2022' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1664912077112 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 30 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.02283E+00  1.02800E+00  1.02219E+00  1.02498E+00  1.04867E+00  1.05263E+00  1.02384E+00  7.91773E-01  7.91122E-01  1.02798E+00  1.05287E+00  1.02005E+00  7.89558E-01  1.02426E+00  1.05196E+00  1.03010E+00  1.01629E+00  1.04779E+00  1.02655E+00  1.05345E+00  7.90504E-01  1.05409E+00  1.02574E+00  1.01319E+00  1.01636E+00  1.02654E+00  1.02312E+00  1.02609E+00  1.02508E+00  1.05239E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 35])  = '/opt/xsdata/xsdir_files/xsdata_endf' ;
DECAY_DATA_FILE_PATH      (idx, [1: 48])  = '/opt/xsdata/JEFF33/Decay_Data/JEFF33-rdd_all.asc' ;
SFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
NFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 3.7E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.10661E-01 0.00051  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  8.89339E-01 6.4E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.79010E-01 0.00014  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  5.25298E-01 0.00012  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.78928E+00 0.00027  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  1.68602E+01 0.00035  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  1.60739E+01 0.00037  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.45258E+01 0.00051  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.36088E+00 0.00035  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 9997390 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.49683E+01 ;
RUNNING_TIME              (idx, 1)        =  6.65283E-01 ;
INIT_TIME                 (idx, [1:  2])  = [  4.35167E-02  4.35167E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.68333E-03  1.68333E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  6.19433E-01  6.19433E-01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  6.64300E-01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 22.49919 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.41043E+01 0.00057 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  7.00829E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128701.35 ;
ALLOC_MEMSIZE             (idx, 1)        = 1400.32;
MEMSIZE                   (idx, 1)        = 1091.76;
XS_MEMSIZE                (idx, 1)        = 387.35;
MAT_MEMSIZE               (idx, 1)        = 27.72;
RES_MEMSIZE               (idx, 1)        = 8.37;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 668.32;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 308.56;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 23 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 241498 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 2 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 14 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 14 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 309 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.99831E-06 0.00013  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.62158E-01 0.00173 ];
U235_FISS                 (idx, [1:   4]) = [  9.00204E-02 0.00098  9.12290E-01 0.00029 ];
U238_FISS                 (idx, [1:   4]) = [  8.65552E-03 0.00345  8.77097E-02 0.00305 ];
U235_CAPT                 (idx, [1:   4]) = [  2.86726E-02 0.00179  2.49119E-01 0.00155 ];
U238_CAPT                 (idx, [1:   4]) = [  4.29957E-02 0.00150  3.73560E-01 0.00112 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 9997390 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 3.65920E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 1150434 1.15117E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 986249 9.86921E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 7860707 7.86557E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -8.17701E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  3.20185E-12 0.00076 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  2.44351E-01 0.00075 ];
TOT_FISSRATE              (idx, [1:   2]) = [  9.86052E-02 0.00076 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  1.14971E-01 0.00061 ];
TOT_ABSRATE               (idx, [1:   2]) = [  2.13577E-01 0.00060 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.99831E-01 0.00013 ];
TOT_FLUX                  (idx, [1:   2]) = [  3.83988E+01 0.00029 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.86423E-01 0.00016 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  1.60777E+01 0.00034 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.86402E+00 0.00073 ];
SIX_FF_F                  (idx, [1:   2]) = [  8.52868E-01 0.00049 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.82670E-01 0.00110 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.55022E+00 0.00120 ];
SIX_FF_LF                 (idx, [1:   2]) = [  2.80565E-01 0.00059 ];
SIX_FF_LT                 (idx, [1:   2]) = [  7.60762E-01 0.00034 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14585E+00 0.00080 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  2.44572E-01 0.00097 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.47807E+00 1.6E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02671E+02 1.7E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  2.44561E-01 0.00098  2.42443E-01 0.00097  2.12865E-03 0.01050 ];
IMP_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
COL_KEFF                  (idx, [1:   2]) = [  2.44393E-01 0.00080 ];
ABS_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14599E+00 0.00034 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.31676E+01 0.00048 ];
IMP_ALF                   (idx, [1:   2]) = [  1.31807E+01 0.00027 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  3.83063E-05 0.00631 ];
IMP_EALF                  (idx, [1:   2]) = [  3.77563E-05 0.00355 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  3.83173E-01 0.00272 ];
IMP_AFGE                  (idx, [1:   2]) = [  3.82479E-01 0.00091 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  3.01927E-02 0.00419  9.24995E-04 0.02176  5.03528E-03 0.00845  4.99499E-03 0.00923  1.16770E-02 0.00547  5.34205E-03 0.00867  2.21837E-03 0.01298 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  5.11675E-01 0.00488  1.33611E-02 0.00014  3.25403E-02 0.00013  1.21199E-01 7.0E-05  3.07058E-01 0.00018  8.66487E-01 0.00026  2.90811E+00 0.00043 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  8.81649E-03 0.00816  2.72268E-04 0.05149  1.50312E-03 0.01995  1.45803E-03 0.02070  3.35929E-03 0.01468  1.58926E-03 0.01893  6.34525E-04 0.03381 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  5.08559E-01 0.01310  1.33544E-02 0.00022  3.25614E-02 0.00032  1.21174E-01 0.00017  3.06933E-01 0.00041  8.66936E-01 0.00062  2.90992E+00 0.00097 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  4.85774E-05 0.00236  4.85134E-05 0.00235  5.61001E-05 0.02312 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.18788E-05 0.00207  1.18632E-05 0.00207  1.37153E-05 0.02293 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  8.68974E-03 0.01043  2.70582E-04 0.06334  1.47407E-03 0.02434  1.47764E-03 0.02685  3.28101E-03 0.01582  1.55596E-03 0.02523  6.30468E-04 0.04396 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  5.07907E-01 0.01635  1.33534E-02 0.00034  3.25553E-02 0.00040  1.21201E-01 0.00021  3.07118E-01 0.00056  8.67255E-01 0.00084  2.90641E+00 0.00138 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.85141E-05 0.00658  4.84099E-05 0.00652  5.75531E-05 0.05433 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.18632E-05 0.00647  1.18378E-05 0.00641  1.40664E-05 0.05415 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  8.66956E-03 0.04306  4.20749E-04 0.22546  1.41031E-03 0.09430  1.20524E-03 0.09711  3.31507E-03 0.07513  1.38061E-03 0.11194  9.37594E-04 0.14050 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  5.79172E-01 0.06119  1.33570E-02 0.00109  3.25956E-02 0.00116  1.21224E-01 0.00076  3.06904E-01 0.00173  8.66676E-01 0.00283  2.90741E+00 0.00357 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  8.72988E-03 0.04178  4.15360E-04 0.22586  1.38787E-03 0.09259  1.24556E-03 0.09633  3.32720E-03 0.07154  1.40820E-03 0.11097  9.45694E-04 0.13655 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  5.86308E-01 0.06010  1.33570E-02 0.00109  3.25985E-02 0.00115  1.21217E-01 0.00074  3.06821E-01 0.00170  8.66841E-01 0.00284  2.90677E+00 0.00356 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.80255E+02 0.04373 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.84821E-05 0.00148 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.18556E-05 0.00104 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  8.55066E-03 0.00690 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.76396E+02 0.00696 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  8.86765E-08 0.00092 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.72310E-06 0.00091  3.72311E-06 0.00093  3.71805E-06 0.00789 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.86564E-05 0.00093  1.86590E-05 0.00094  1.84121E-05 0.00995 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  1.28519E-01 0.00095  1.31202E-01 0.00095  4.24282E-02 0.00932 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.01810E+01 0.00891 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  1.60739E+01 0.00037  2.57880E+01 0.00064 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '1' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.19768E+05 0.00592  5.88059E+05 0.00288  1.34958E+06 0.00106  1.86066E+06 0.00105  1.84977E+06 0.00136  1.85489E+06 0.00174  1.11670E+06 0.00071  9.57746E+05 0.00095  1.10399E+06 0.00058  8.77001E+05 0.00170  7.14629E+05 0.00087  5.83715E+05 0.00143  5.16665E+05 0.00170  4.27428E+05 0.00254  3.92733E+05 0.00128  2.75326E+05 0.00171  1.39701E+05 0.00133  4.28947E+05 0.00167  3.47953E+05 0.00253  5.94650E+05 0.00138  5.32748E+05 0.00166  3.49475E+05 0.00090  1.92693E+05 0.00325  2.02975E+05 0.00278  1.76356E+05 0.00268  1.46303E+05 0.00098  2.24358E+05 0.00450  4.83019E+04 0.00619  5.92105E+04 0.00450  5.26946E+04 0.00815  2.94514E+04 0.00620  4.97290E+04 0.00379  3.25956E+04 0.00419  2.62304E+04 0.00305  4.86707E+03 0.01040  4.75796E+03 0.01172  4.88607E+03 0.00679  4.92303E+03 0.00236  4.89904E+03 0.01642  4.83147E+03 0.01171  4.94243E+03 0.01241  4.70643E+03 0.01611  8.65710E+03 0.00609  1.37691E+04 0.00665  1.72633E+04 0.00535  4.52280E+04 0.00390  4.71173E+04 0.00323  5.08137E+04 0.00385  3.28915E+04 0.00579  2.31289E+04 0.00643  1.72313E+04 0.00472  1.90120E+04 0.00445  3.22346E+04 0.00382  3.86646E+04 0.00242  6.20163E+04 0.00296  7.52429E+04 0.00222  8.62588E+04 0.00365  4.48389E+04 0.00372  2.84710E+04 0.00226  1.89317E+04 0.00196  1.59734E+04 0.00479  1.47565E+04 0.00580  1.17467E+04 0.00310  7.66209E+03 0.00661  6.84056E+03 0.00506  5.88445E+03 0.00390  4.82432E+03 0.00292  3.62930E+03 0.00784  2.31683E+03 0.00840  7.98255E+02 0.00493 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.00262E+00 0.00110 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  9.17183E+00 0.00081  3.25588E-01 0.00202 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.14060E-01 8.5E-05  7.48975E-01 0.00019 ];
INF_CAPT                  (idx, [1:   4]) = [  3.34802E-03 0.00063  1.82743E-02 0.00252 ];
INF_ABS                   (idx, [1:   4]) = [  4.93345E-03 0.00055  5.01808E-02 0.00192 ];
INF_FISS                  (idx, [1:   4]) = [  1.58543E-03 0.00086  3.19064E-02 0.00165 ];
INF_NSF                   (idx, [1:   4]) = [  3.96966E-03 0.00085  7.77464E-02 0.00165 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.50383E+00 4.2E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02921E+02 3.6E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  2.85387E-08 0.00182  1.97765E-06 0.00023 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  4.09128E-01 9.3E-05  6.98581E-01 0.00029 ];
INF_SCATT1                (idx, [1:   4]) = [  1.21004E-01 0.00040  2.06769E-01 0.00064 ];
INF_SCATT2                (idx, [1:   4]) = [  4.75618E-02 0.00067  7.65010E-02 0.00318 ];
INF_SCATT3                (idx, [1:   4]) = [  4.18367E-03 0.00518  3.05834E-02 0.00573 ];
INF_SCATT4                (idx, [1:   4]) = [ -4.01901E-03 0.00612  1.40776E-02 0.01399 ];
INF_SCATT5                (idx, [1:   4]) = [  9.05025E-05 0.22242  7.40571E-03 0.01253 ];
INF_SCATT6                (idx, [1:   4]) = [  2.23196E-03 0.00663  4.46781E-03 0.01880 ];
INF_SCATT7                (idx, [1:   4]) = [  3.47576E-04 0.04156  3.11727E-03 0.03022 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  4.09137E-01 9.3E-05  6.98581E-01 0.00029 ];
INF_SCATTP1               (idx, [1:   4]) = [  1.21004E-01 0.00040  2.06769E-01 0.00064 ];
INF_SCATTP2               (idx, [1:   4]) = [  4.75618E-02 0.00067  7.65010E-02 0.00318 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.18379E-03 0.00517  3.05834E-02 0.00573 ];
INF_SCATTP4               (idx, [1:   4]) = [ -4.01873E-03 0.00612  1.40776E-02 0.01399 ];
INF_SCATTP5               (idx, [1:   4]) = [  9.05274E-05 0.22213  7.40571E-03 0.01253 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.23191E-03 0.00662  4.46781E-03 0.01880 ];
INF_SCATTP7               (idx, [1:   4]) = [  3.47461E-04 0.04159  3.11727E-03 0.03022 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.36534E-01 0.00025  5.05778E-01 0.00056 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.40924E+00 0.00025  6.59052E-01 0.00056 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  4.92432E-03 0.00053  5.01808E-02 0.00192 ];
INF_REMXS                 (idx, [1:   4]) = [  8.74977E-03 0.00167  5.50020E-02 0.00168 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37701E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49395E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.94297E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.92615E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55062E-04 ;
I135_BR                   (idx, 1)        =  8.34918E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  4.05311E-01 0.00011  3.81755E-03 0.00170  4.60735E-03 0.00656  6.93973E-01 0.00031 ];
INF_S1                    (idx, [1:   8]) = [  1.19856E-01 0.00039  1.14820E-03 0.00236  1.60596E-03 0.01397  2.05164E-01 0.00074 ];
INF_S2                    (idx, [1:   8]) = [  4.78208E-02 0.00066 -2.58992E-04 0.00755  7.88430E-04 0.03426  7.57125E-02 0.00308 ];
INF_S3                    (idx, [1:   8]) = [  4.59635E-03 0.00440 -4.12674E-04 0.00619  2.56405E-04 0.06385  3.03270E-02 0.00555 ];
INF_S4                    (idx, [1:   8]) = [ -3.86181E-03 0.00608 -1.57205E-04 0.00959  9.51140E-06 1.00000  1.40681E-02 0.01413 ];
INF_S5                    (idx, [1:   8]) = [  1.13053E-04 0.17403 -2.25506E-05 0.04127 -1.01666E-04 0.23171  7.50738E-03 0.01420 ];
INF_S6                    (idx, [1:   8]) = [  2.23789E-03 0.00686 -5.92835E-06 0.23797 -9.00526E-05 0.15882  4.55786E-03 0.01735 ];
INF_S7                    (idx, [1:   8]) = [  3.56696E-04 0.04124 -9.12009E-06 0.07641 -8.08998E-05 0.13903  3.19817E-03 0.02967 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  4.05320E-01 0.00011  3.81755E-03 0.00170  4.60735E-03 0.00656  6.93973E-01 0.00031 ];
INF_SP1                   (idx, [1:   8]) = [  1.19856E-01 0.00039  1.14820E-03 0.00236  1.60596E-03 0.01397  2.05164E-01 0.00074 ];
INF_SP2                   (idx, [1:   8]) = [  4.78208E-02 0.00066 -2.58992E-04 0.00755  7.88430E-04 0.03426  7.57125E-02 0.00308 ];
INF_SP3                   (idx, [1:   8]) = [  4.59647E-03 0.00440 -4.12674E-04 0.00619  2.56405E-04 0.06385  3.03270E-02 0.00555 ];
INF_SP4                   (idx, [1:   8]) = [ -3.86153E-03 0.00609 -1.57205E-04 0.00959  9.51140E-06 1.00000  1.40681E-02 0.01413 ];
INF_SP5                   (idx, [1:   8]) = [  1.13078E-04 0.17380 -2.25506E-05 0.04127 -1.01666E-04 0.23171  7.50738E-03 0.01420 ];
INF_SP6                   (idx, [1:   8]) = [  2.23784E-03 0.00685 -5.92835E-06 0.23797 -9.00526E-05 0.15882  4.55786E-03 0.01735 ];
INF_SP7                   (idx, [1:   8]) = [  3.56581E-04 0.04127 -9.12009E-06 0.07641 -8.08998E-05 0.13903  3.19817E-03 0.02967 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  4.71516E-01 0.00142  1.54405E+00 0.00920 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  6.67600E-01 0.00232  1.94402E+00 0.01773 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  6.72540E-01 0.00255 -1.95641E+00 0.03028 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.96068E-01 0.00148  5.15355E-01 0.01233 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  7.06945E-01 0.00142  2.15955E-01 0.00915 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  4.99312E-01 0.00232  1.71676E-01 0.01726 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  4.95646E-01 0.00254 -1.71018E-01 0.03083 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.12588E+00 0.00149  6.47208E-01 0.01266 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  8.74014E-03 0.01641  2.69708E-04 0.09032  1.53583E-03 0.04396  1.41331E-03 0.03760  3.39130E-03 0.02740  1.53077E-03 0.03334  5.99222E-04 0.06101 ];
LAMBDA                    (idx, [1:  14]) = [  4.99440E-01 0.02307  1.33541E-02 0.00035  3.25485E-02 0.00054  1.21148E-01 0.00030  3.06936E-01 0.00075  8.66517E-01 0.00112  2.90671E+00 0.00172 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'May 19 2021 21:58:01' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 27])  = 'hexassembly_2D_absorber.inp' ;
WORKING_DIRECTORY         (idx, [1: 26])  = '/nfs/home/jgj/projects/UPV' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode01' ;
CPU_TYPE                  (idx, [1: 46])  = 'AMD Ryzen Threadripper 3970X 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 137367609.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct  4 21:34:37 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct  4 21:35:17 2022' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1664912077112 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 30 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.02283E+00  1.02800E+00  1.02219E+00  1.02498E+00  1.04867E+00  1.05263E+00  1.02384E+00  7.91773E-01  7.91122E-01  1.02798E+00  1.05287E+00  1.02005E+00  7.89558E-01  1.02426E+00  1.05196E+00  1.03010E+00  1.01629E+00  1.04779E+00  1.02655E+00  1.05345E+00  7.90504E-01  1.05409E+00  1.02574E+00  1.01319E+00  1.01636E+00  1.02654E+00  1.02312E+00  1.02609E+00  1.02508E+00  1.05239E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 35])  = '/opt/xsdata/xsdir_files/xsdata_endf' ;
DECAY_DATA_FILE_PATH      (idx, [1: 48])  = '/opt/xsdata/JEFF33/Decay_Data/JEFF33-rdd_all.asc' ;
SFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
NFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 3.7E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.10661E-01 0.00051  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  8.89339E-01 6.4E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.79010E-01 0.00014  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  5.25298E-01 0.00012  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.78928E+00 0.00027  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  1.68602E+01 0.00035  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  1.60739E+01 0.00037  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.45258E+01 0.00051  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.36088E+00 0.00035  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 9997390 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.49683E+01 ;
RUNNING_TIME              (idx, 1)        =  6.65300E-01 ;
INIT_TIME                 (idx, [1:  2])  = [  4.35167E-02  4.35167E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.68333E-03  1.68333E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  6.19433E-01  6.19433E-01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  6.64300E-01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 22.49862 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.41043E+01 0.00057 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  7.00812E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128701.35 ;
ALLOC_MEMSIZE             (idx, 1)        = 1400.32;
MEMSIZE                   (idx, 1)        = 1091.76;
XS_MEMSIZE                (idx, 1)        = 387.35;
MAT_MEMSIZE               (idx, 1)        = 27.72;
RES_MEMSIZE               (idx, 1)        = 8.37;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 668.32;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 308.56;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 23 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 241498 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 2 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 14 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 14 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 309 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.99831E-06 0.00013  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.62158E-01 0.00173 ];
U235_FISS                 (idx, [1:   4]) = [  9.00204E-02 0.00098  9.12290E-01 0.00029 ];
U238_FISS                 (idx, [1:   4]) = [  8.65552E-03 0.00345  8.77097E-02 0.00305 ];
U235_CAPT                 (idx, [1:   4]) = [  2.86726E-02 0.00179  2.49119E-01 0.00155 ];
U238_CAPT                 (idx, [1:   4]) = [  4.29957E-02 0.00150  3.73560E-01 0.00112 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 9997390 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 3.65920E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 1150434 1.15117E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 986249 9.86921E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 7860707 7.86557E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -8.17701E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  3.20185E-12 0.00076 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  2.44351E-01 0.00075 ];
TOT_FISSRATE              (idx, [1:   2]) = [  9.86052E-02 0.00076 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  1.14971E-01 0.00061 ];
TOT_ABSRATE               (idx, [1:   2]) = [  2.13577E-01 0.00060 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.99831E-01 0.00013 ];
TOT_FLUX                  (idx, [1:   2]) = [  3.83988E+01 0.00029 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.86423E-01 0.00016 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  1.60777E+01 0.00034 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.86402E+00 0.00073 ];
SIX_FF_F                  (idx, [1:   2]) = [  8.52868E-01 0.00049 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.82670E-01 0.00110 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.55022E+00 0.00120 ];
SIX_FF_LF                 (idx, [1:   2]) = [  2.80565E-01 0.00059 ];
SIX_FF_LT                 (idx, [1:   2]) = [  7.60762E-01 0.00034 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14585E+00 0.00080 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  2.44572E-01 0.00097 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.47807E+00 1.6E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02671E+02 1.7E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  2.44561E-01 0.00098  2.42443E-01 0.00097  2.12865E-03 0.01050 ];
IMP_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
COL_KEFF                  (idx, [1:   2]) = [  2.44393E-01 0.00080 ];
ABS_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14599E+00 0.00034 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.31676E+01 0.00048 ];
IMP_ALF                   (idx, [1:   2]) = [  1.31807E+01 0.00027 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  3.83063E-05 0.00631 ];
IMP_EALF                  (idx, [1:   2]) = [  3.77563E-05 0.00355 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  3.83173E-01 0.00272 ];
IMP_AFGE                  (idx, [1:   2]) = [  3.82479E-01 0.00091 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  3.01927E-02 0.00419  9.24995E-04 0.02176  5.03528E-03 0.00845  4.99499E-03 0.00923  1.16770E-02 0.00547  5.34205E-03 0.00867  2.21837E-03 0.01298 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  5.11675E-01 0.00488  1.33611E-02 0.00014  3.25403E-02 0.00013  1.21199E-01 7.0E-05  3.07058E-01 0.00018  8.66487E-01 0.00026  2.90811E+00 0.00043 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  8.81649E-03 0.00816  2.72268E-04 0.05149  1.50312E-03 0.01995  1.45803E-03 0.02070  3.35929E-03 0.01468  1.58926E-03 0.01893  6.34525E-04 0.03381 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  5.08559E-01 0.01310  1.33544E-02 0.00022  3.25614E-02 0.00032  1.21174E-01 0.00017  3.06933E-01 0.00041  8.66936E-01 0.00062  2.90992E+00 0.00097 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  4.85774E-05 0.00236  4.85134E-05 0.00235  5.61001E-05 0.02312 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.18788E-05 0.00207  1.18632E-05 0.00207  1.37153E-05 0.02293 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  8.68974E-03 0.01043  2.70582E-04 0.06334  1.47407E-03 0.02434  1.47764E-03 0.02685  3.28101E-03 0.01582  1.55596E-03 0.02523  6.30468E-04 0.04396 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  5.07907E-01 0.01635  1.33534E-02 0.00034  3.25553E-02 0.00040  1.21201E-01 0.00021  3.07118E-01 0.00056  8.67255E-01 0.00084  2.90641E+00 0.00138 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.85141E-05 0.00658  4.84099E-05 0.00652  5.75531E-05 0.05433 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.18632E-05 0.00647  1.18378E-05 0.00641  1.40664E-05 0.05415 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  8.66956E-03 0.04306  4.20749E-04 0.22546  1.41031E-03 0.09430  1.20524E-03 0.09711  3.31507E-03 0.07513  1.38061E-03 0.11194  9.37594E-04 0.14050 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  5.79172E-01 0.06119  1.33570E-02 0.00109  3.25956E-02 0.00116  1.21224E-01 0.00076  3.06904E-01 0.00173  8.66676E-01 0.00283  2.90741E+00 0.00357 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  8.72988E-03 0.04178  4.15360E-04 0.22586  1.38787E-03 0.09259  1.24556E-03 0.09633  3.32720E-03 0.07154  1.40820E-03 0.11097  9.45694E-04 0.13655 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  5.86308E-01 0.06010  1.33570E-02 0.00109  3.25985E-02 0.00115  1.21217E-01 0.00074  3.06821E-01 0.00170  8.66841E-01 0.00284  2.90677E+00 0.00356 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.80255E+02 0.04373 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.84821E-05 0.00148 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.18556E-05 0.00104 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  8.55066E-03 0.00690 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.76396E+02 0.00696 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  8.86765E-08 0.00092 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.72310E-06 0.00091  3.72311E-06 0.00093  3.71805E-06 0.00789 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.86564E-05 0.00093  1.86590E-05 0.00094  1.84121E-05 0.00995 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  1.28519E-01 0.00095  1.31202E-01 0.00095  4.24282E-02 0.00932 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.01810E+01 0.00891 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  1.60739E+01 0.00037  2.57880E+01 0.00064 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '2' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  6.36997E+04 0.00479  3.04364E+05 0.00152  7.01161E+05 0.00099  9.61194E+05 0.00218  9.42273E+05 0.00249  9.35548E+05 0.00091  5.58805E+05 0.00099  4.76324E+05 0.00142  5.56826E+05 0.00223  4.42059E+05 0.00068  3.61627E+05 0.00255  2.94842E+05 0.00346  2.59749E+05 0.00312  2.12314E+05 0.00184  1.95055E+05 0.00241  1.35778E+05 0.00239  6.79079E+04 0.00232  2.09873E+05 0.00214  1.72484E+05 0.00217  2.91638E+05 0.00254  2.56329E+05 0.00415  1.67541E+05 0.00239  9.14500E+04 0.00242  9.49635E+04 0.00367  8.21140E+04 0.00118  6.80339E+04 0.00368  1.02733E+05 0.00158  2.22195E+04 0.00555  2.68685E+04 0.00631  2.41182E+04 0.00626  1.34166E+04 0.00797  2.26498E+04 0.00442  1.50367E+04 0.00328  1.20526E+04 0.00519  2.17189E+03 0.01156  2.12733E+03 0.00492  2.19838E+03 0.02073  2.25268E+03 0.01690  2.17132E+03 0.01100  2.20251E+03 0.01172  2.25095E+03 0.01408  2.00990E+03 0.01317  3.94664E+03 0.01800  6.34246E+03 0.00726  7.94561E+03 0.01171  2.01633E+04 0.00388  2.08067E+04 0.00554  2.21261E+04 0.00725  1.42170E+04 0.00464  9.79612E+03 0.00560  7.24563E+03 0.00710  8.00089E+03 0.00913  1.38251E+04 0.00621  1.61143E+04 0.00654  2.58411E+04 0.00250  3.12862E+04 0.00396  3.56901E+04 0.00540  1.84185E+04 0.00687  1.17828E+04 0.00668  7.74647E+03 0.00915  6.57191E+03 0.00923  6.12085E+03 0.00608  4.85986E+03 0.00398  3.08338E+03 0.00418  2.76962E+03 0.01932  2.41888E+03 0.00672  2.02320E+03 0.01091  1.48454E+03 0.01326  9.53335E+02 0.00954  3.26186E+02 0.01922 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.19691E+00 0.00127 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.59864E+00 0.00127  1.36731E-01 0.00340 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.06806E-01 0.00024  7.41144E-01 0.00067 ];
INF_CAPT                  (idx, [1:   4]) = [  2.25091E-03 0.00060  1.69138E-02 0.00116 ];
INF_ABS                   (idx, [1:   4]) = [  3.83443E-03 0.00043  4.98098E-02 0.00220 ];
INF_FISS                  (idx, [1:   4]) = [  1.58352E-03 0.00058  3.28961E-02 0.00276 ];
INF_NSF                   (idx, [1:   4]) = [  3.96798E-03 0.00053  8.01579E-02 0.00276 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.50580E+00 5.1E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02939E+02 3.8E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  2.64627E-08 0.00045  1.96265E-06 0.00053 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  4.02967E-01 0.00023  6.91456E-01 0.00087 ];
INF_SCATT1                (idx, [1:   4]) = [  1.17810E-01 0.00046  2.03513E-01 0.00113 ];
INF_SCATT2                (idx, [1:   4]) = [  4.64779E-02 0.00115  7.48421E-02 0.00353 ];
INF_SCATT3                (idx, [1:   4]) = [  4.14934E-03 0.01197  2.92806E-02 0.00779 ];
INF_SCATT4                (idx, [1:   4]) = [ -3.86163E-03 0.02016  1.36522E-02 0.01366 ];
INF_SCATT5                (idx, [1:   4]) = [  4.72543E-05 0.66614  6.91319E-03 0.02621 ];
INF_SCATT6                (idx, [1:   4]) = [  2.12159E-03 0.01009  3.95524E-03 0.05150 ];
INF_SCATT7                (idx, [1:   4]) = [  3.19101E-04 0.02474  2.96831E-03 0.06705 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  4.02977E-01 0.00023  6.91456E-01 0.00087 ];
INF_SCATTP1               (idx, [1:   4]) = [  1.17810E-01 0.00046  2.03513E-01 0.00113 ];
INF_SCATTP2               (idx, [1:   4]) = [  4.64779E-02 0.00115  7.48421E-02 0.00353 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.14935E-03 0.01194  2.92806E-02 0.00779 ];
INF_SCATTP4               (idx, [1:   4]) = [ -3.86139E-03 0.02017  1.36522E-02 0.01366 ];
INF_SCATTP5               (idx, [1:   4]) = [  4.73523E-05 0.66726  6.91319E-03 0.02621 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.12151E-03 0.01006  3.95524E-03 0.05150 ];
INF_SCATTP7               (idx, [1:   4]) = [  3.19169E-04 0.02477  2.96831E-03 0.06705 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.32938E-01 0.00029  5.01321E-01 0.00082 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.43099E+00 0.00029  6.64912E-01 0.00082 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  3.82436E-03 0.00036  4.98098E-02 0.00220 ];
INF_REMXS                 (idx, [1:   4]) = [  7.22159E-03 0.00191  5.43826E-02 0.00333 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37701E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49395E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.94297E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.92615E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55062E-04 ;
I135_BR                   (idx, 1)        =  8.34918E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  3.99585E-01 0.00022  3.38246E-03 0.00213  4.69462E-03 0.02317  6.86761E-01 0.00095 ];
INF_S1                    (idx, [1:   8]) = [  1.16793E-01 0.00043  1.01701E-03 0.00555  1.65540E-03 0.01855  2.01858E-01 0.00114 ];
INF_S2                    (idx, [1:   8]) = [  4.67053E-02 0.00116 -2.27352E-04 0.00991  8.31597E-04 0.04099  7.40105E-02 0.00340 ];
INF_S3                    (idx, [1:   8]) = [  4.51162E-03 0.01091 -3.62280E-04 0.00660  2.85303E-04 0.13454  2.89953E-02 0.00715 ];
INF_S4                    (idx, [1:   8]) = [ -3.72260E-03 0.02057 -1.39026E-04 0.01040  1.60582E-05 1.00000  1.36361E-02 0.01482 ];
INF_S5                    (idx, [1:   8]) = [  7.07294E-05 0.42311 -2.34751E-05 0.11202 -6.89857E-05 0.31536  6.98218E-03 0.02592 ];
INF_S6                    (idx, [1:   8]) = [  2.12863E-03 0.00963 -7.04395E-06 0.20024 -8.21000E-05 0.29118  4.03734E-03 0.04851 ];
INF_S7                    (idx, [1:   8]) = [  3.27499E-04 0.01770 -8.39829E-06 0.30820 -5.45535E-05 0.36473  3.02287E-03 0.06197 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.99595E-01 0.00022  3.38246E-03 0.00213  4.69462E-03 0.02317  6.86761E-01 0.00095 ];
INF_SP1                   (idx, [1:   8]) = [  1.16793E-01 0.00044  1.01701E-03 0.00555  1.65540E-03 0.01855  2.01858E-01 0.00114 ];
INF_SP2                   (idx, [1:   8]) = [  4.67052E-02 0.00116 -2.27352E-04 0.00991  8.31597E-04 0.04099  7.40105E-02 0.00340 ];
INF_SP3                   (idx, [1:   8]) = [  4.51163E-03 0.01088 -3.62280E-04 0.00660  2.85303E-04 0.13454  2.89953E-02 0.00715 ];
INF_SP4                   (idx, [1:   8]) = [ -3.72236E-03 0.02058 -1.39026E-04 0.01040  1.60582E-05 1.00000  1.36361E-02 0.01482 ];
INF_SP5                   (idx, [1:   8]) = [  7.08274E-05 0.42424 -2.34751E-05 0.11202 -6.89857E-05 0.31536  6.98218E-03 0.02592 ];
INF_SP6                   (idx, [1:   8]) = [  2.12855E-03 0.00960 -7.04395E-06 0.20024 -8.21000E-05 0.29118  4.03734E-03 0.04851 ];
INF_SP7                   (idx, [1:   8]) = [  3.27567E-04 0.01772 -8.39829E-06 0.30820 -5.45535E-05 0.36473  3.02287E-03 0.06197 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  1.37388E-01 0.00119  3.34656E-01 0.01019 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  1.56044E-01 0.00129  4.38260E-01 0.00999 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  1.94531E-01 0.00184  1.04405E+01 0.89581 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  9.72115E-02 0.00195  1.52202E-01 0.01275 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  2.42624E+00 0.00119  9.96457E-01 0.01004 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  2.13617E+00 0.00129  7.60889E-01 0.01004 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.71354E+00 0.00183  3.69803E-02 0.62965 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  3.42900E+00 0.00195  2.19150E+00 0.01283 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  8.76628E-03 0.02629  2.78555E-04 0.13161  1.43159E-03 0.05707  1.39223E-03 0.05524  3.35291E-03 0.04174  1.68331E-03 0.05816  6.27692E-04 0.08806 ];
LAMBDA                    (idx, [1:  14]) = [  5.16411E-01 0.03397  1.33638E-02 0.00054  3.25672E-02 0.00068  1.21239E-01 0.00044  3.06609E-01 0.00095  8.66219E-01 0.00145  2.91721E+00 0.00202 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'May 19 2021 21:58:01' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 27])  = 'hexassembly_2D_absorber.inp' ;
WORKING_DIRECTORY         (idx, [1: 26])  = '/nfs/home/jgj/projects/UPV' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode01' ;
CPU_TYPE                  (idx, [1: 46])  = 'AMD Ryzen Threadripper 3970X 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 137367609.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct  4 21:34:37 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct  4 21:35:17 2022' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1664912077112 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 30 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.02283E+00  1.02800E+00  1.02219E+00  1.02498E+00  1.04867E+00  1.05263E+00  1.02384E+00  7.91773E-01  7.91122E-01  1.02798E+00  1.05287E+00  1.02005E+00  7.89558E-01  1.02426E+00  1.05196E+00  1.03010E+00  1.01629E+00  1.04779E+00  1.02655E+00  1.05345E+00  7.90504E-01  1.05409E+00  1.02574E+00  1.01319E+00  1.01636E+00  1.02654E+00  1.02312E+00  1.02609E+00  1.02508E+00  1.05239E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 35])  = '/opt/xsdata/xsdir_files/xsdata_endf' ;
DECAY_DATA_FILE_PATH      (idx, [1: 48])  = '/opt/xsdata/JEFF33/Decay_Data/JEFF33-rdd_all.asc' ;
SFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
NFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 3.7E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.10661E-01 0.00051  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  8.89339E-01 6.4E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.79010E-01 0.00014  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  5.25298E-01 0.00012  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.78928E+00 0.00027  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  1.68602E+01 0.00035  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  1.60739E+01 0.00037  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.45258E+01 0.00051  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.36088E+00 0.00035  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 9997390 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.49683E+01 ;
RUNNING_TIME              (idx, 1)        =  6.65300E-01 ;
INIT_TIME                 (idx, [1:  2])  = [  4.35167E-02  4.35167E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.68333E-03  1.68333E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  6.19433E-01  6.19433E-01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  6.64300E-01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 22.49862 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.41043E+01 0.00057 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  7.00812E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128701.35 ;
ALLOC_MEMSIZE             (idx, 1)        = 1400.32;
MEMSIZE                   (idx, 1)        = 1091.76;
XS_MEMSIZE                (idx, 1)        = 387.35;
MAT_MEMSIZE               (idx, 1)        = 27.72;
RES_MEMSIZE               (idx, 1)        = 8.37;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 668.32;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 308.56;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 23 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 241498 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 2 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 14 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 14 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 309 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.99831E-06 0.00013  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.62158E-01 0.00173 ];
U235_FISS                 (idx, [1:   4]) = [  9.00204E-02 0.00098  9.12290E-01 0.00029 ];
U238_FISS                 (idx, [1:   4]) = [  8.65552E-03 0.00345  8.77097E-02 0.00305 ];
U235_CAPT                 (idx, [1:   4]) = [  2.86726E-02 0.00179  2.49119E-01 0.00155 ];
U238_CAPT                 (idx, [1:   4]) = [  4.29957E-02 0.00150  3.73560E-01 0.00112 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 9997390 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 3.65920E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 1150434 1.15117E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 986249 9.86921E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 7860707 7.86557E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -8.17701E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  3.20185E-12 0.00076 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  2.44351E-01 0.00075 ];
TOT_FISSRATE              (idx, [1:   2]) = [  9.86052E-02 0.00076 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  1.14971E-01 0.00061 ];
TOT_ABSRATE               (idx, [1:   2]) = [  2.13577E-01 0.00060 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.99831E-01 0.00013 ];
TOT_FLUX                  (idx, [1:   2]) = [  3.83988E+01 0.00029 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.86423E-01 0.00016 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  1.60777E+01 0.00034 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.86402E+00 0.00073 ];
SIX_FF_F                  (idx, [1:   2]) = [  8.52868E-01 0.00049 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.82670E-01 0.00110 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.55022E+00 0.00120 ];
SIX_FF_LF                 (idx, [1:   2]) = [  2.80565E-01 0.00059 ];
SIX_FF_LT                 (idx, [1:   2]) = [  7.60762E-01 0.00034 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14585E+00 0.00080 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  2.44572E-01 0.00097 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.47807E+00 1.6E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02671E+02 1.7E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  2.44561E-01 0.00098  2.42443E-01 0.00097  2.12865E-03 0.01050 ];
IMP_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
COL_KEFF                  (idx, [1:   2]) = [  2.44393E-01 0.00080 ];
ABS_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14599E+00 0.00034 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.31676E+01 0.00048 ];
IMP_ALF                   (idx, [1:   2]) = [  1.31807E+01 0.00027 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  3.83063E-05 0.00631 ];
IMP_EALF                  (idx, [1:   2]) = [  3.77563E-05 0.00355 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  3.83173E-01 0.00272 ];
IMP_AFGE                  (idx, [1:   2]) = [  3.82479E-01 0.00091 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  3.01927E-02 0.00419  9.24995E-04 0.02176  5.03528E-03 0.00845  4.99499E-03 0.00923  1.16770E-02 0.00547  5.34205E-03 0.00867  2.21837E-03 0.01298 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  5.11675E-01 0.00488  1.33611E-02 0.00014  3.25403E-02 0.00013  1.21199E-01 7.0E-05  3.07058E-01 0.00018  8.66487E-01 0.00026  2.90811E+00 0.00043 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  8.81649E-03 0.00816  2.72268E-04 0.05149  1.50312E-03 0.01995  1.45803E-03 0.02070  3.35929E-03 0.01468  1.58926E-03 0.01893  6.34525E-04 0.03381 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  5.08559E-01 0.01310  1.33544E-02 0.00022  3.25614E-02 0.00032  1.21174E-01 0.00017  3.06933E-01 0.00041  8.66936E-01 0.00062  2.90992E+00 0.00097 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  4.85774E-05 0.00236  4.85134E-05 0.00235  5.61001E-05 0.02312 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.18788E-05 0.00207  1.18632E-05 0.00207  1.37153E-05 0.02293 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  8.68974E-03 0.01043  2.70582E-04 0.06334  1.47407E-03 0.02434  1.47764E-03 0.02685  3.28101E-03 0.01582  1.55596E-03 0.02523  6.30468E-04 0.04396 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  5.07907E-01 0.01635  1.33534E-02 0.00034  3.25553E-02 0.00040  1.21201E-01 0.00021  3.07118E-01 0.00056  8.67255E-01 0.00084  2.90641E+00 0.00138 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.85141E-05 0.00658  4.84099E-05 0.00652  5.75531E-05 0.05433 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.18632E-05 0.00647  1.18378E-05 0.00641  1.40664E-05 0.05415 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  8.66956E-03 0.04306  4.20749E-04 0.22546  1.41031E-03 0.09430  1.20524E-03 0.09711  3.31507E-03 0.07513  1.38061E-03 0.11194  9.37594E-04 0.14050 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  5.79172E-01 0.06119  1.33570E-02 0.00109  3.25956E-02 0.00116  1.21224E-01 0.00076  3.06904E-01 0.00173  8.66676E-01 0.00283  2.90741E+00 0.00357 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  8.72988E-03 0.04178  4.15360E-04 0.22586  1.38787E-03 0.09259  1.24556E-03 0.09633  3.32720E-03 0.07154  1.40820E-03 0.11097  9.45694E-04 0.13655 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  5.86308E-01 0.06010  1.33570E-02 0.00109  3.25985E-02 0.00115  1.21217E-01 0.00074  3.06821E-01 0.00170  8.66841E-01 0.00284  2.90677E+00 0.00356 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.80255E+02 0.04373 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.84821E-05 0.00148 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.18556E-05 0.00104 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  8.55066E-03 0.00690 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.76396E+02 0.00696 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  8.86765E-08 0.00092 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.72310E-06 0.00091  3.72311E-06 0.00093  3.71805E-06 0.00789 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.86564E-05 0.00093  1.86590E-05 0.00094  1.84121E-05 0.00995 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  1.28519E-01 0.00095  1.31202E-01 0.00095  4.24282E-02 0.00932 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.01810E+01 0.00891 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  1.60739E+01 0.00037  2.57880E+01 0.00064 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '3' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  6.66294E+04 0.00556  3.21948E+05 0.00388  7.35801E+05 0.00221  1.00865E+06 0.00187  9.89310E+05 0.00045  9.77222E+05 0.00257  5.80430E+05 0.00380  4.93584E+05 0.00267  5.74726E+05 0.00302  4.54299E+05 0.00278  3.71213E+05 0.00327  3.01093E+05 0.00316  2.64766E+05 0.00266  2.15600E+05 0.00329  1.97274E+05 0.00422  1.37346E+05 0.00316  6.86433E+04 0.00221  2.11786E+05 0.00407  1.73457E+05 0.00439  2.92534E+05 0.00301  2.60794E+05 0.00251  1.70356E+05 0.00248  9.61351E+04 0.00365  1.00817E+05 0.00283  8.79949E+04 0.00533  7.36998E+04 0.00345  1.14101E+05 0.00384  2.48342E+04 0.00557  3.05079E+04 0.00702  2.75841E+04 0.00399  1.54457E+04 0.00332  2.64465E+04 0.00795  1.73703E+04 0.00542  1.38799E+04 0.00930  2.54619E+03 0.00644  2.51997E+03 0.01332  2.49194E+03 0.01680  2.57856E+03 0.00592  2.57882E+03 0.00929  2.56179E+03 0.01485  2.65772E+03 0.01922  2.44361E+03 0.01427  4.58926E+03 0.01227  7.38706E+03 0.00811  9.31657E+03 0.00643  2.40455E+04 0.00333  2.47223E+04 0.00667  2.66938E+04 0.00268  1.72362E+04 0.00673  1.22193E+04 0.00471  9.01725E+03 0.00906  9.84176E+03 0.00682  1.69403E+04 0.00516  2.01556E+04 0.00290  3.25078E+04 0.00339  3.91993E+04 0.00502  4.47284E+04 0.00395  2.32510E+04 0.00865  1.47787E+04 0.00365  9.78323E+03 0.00496  8.25956E+03 0.00695  7.65846E+03 0.00331  6.12897E+03 0.00854  3.98342E+03 0.00530  3.55132E+03 0.00345  3.05977E+03 0.00889  2.50963E+03 0.01054  1.87885E+03 0.01183  1.18702E+03 0.00942  4.27090E+02 0.01041 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.21608E+00 0.00115 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.78019E+00 0.00143  1.69831E-01 0.00275 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.05417E-01 0.00035  7.44354E-01 0.00055 ];
INF_CAPT                  (idx, [1:   4]) = [  2.28256E-03 0.00122  1.68895E-02 0.00054 ];
INF_ABS                   (idx, [1:   4]) = [  3.89436E-03 0.00091  4.94632E-02 0.00132 ];
INF_FISS                  (idx, [1:   4]) = [  1.61181E-03 0.00088  3.25737E-02 0.00175 ];
INF_NSF                   (idx, [1:   4]) = [  4.04016E-03 0.00085  7.93723E-02 0.00175 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.50661E+00 7.5E-05  2.43670E+00 5.9E-09 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02947E+02 9.5E-06  2.02270E+02 5.9E-09 ];
INF_INVV                  (idx, [1:   4]) = [  2.80986E-08 0.00155  1.97442E-06 0.00066 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  4.01519E-01 0.00035  6.95026E-01 0.00055 ];
INF_SCATT1                (idx, [1:   4]) = [  1.17236E-01 0.00057  2.04895E-01 0.00154 ];
INF_SCATT2                (idx, [1:   4]) = [  4.63201E-02 0.00125  7.54482E-02 0.00248 ];
INF_SCATT3                (idx, [1:   4]) = [  4.34427E-03 0.01570  2.98555E-02 0.00631 ];
INF_SCATT4                (idx, [1:   4]) = [ -3.66924E-03 0.01416  1.38752E-02 0.01573 ];
INF_SCATT5                (idx, [1:   4]) = [  1.87338E-04 0.29484  7.38054E-03 0.02039 ];
INF_SCATT6                (idx, [1:   4]) = [  2.10696E-03 0.01149  4.47094E-03 0.02327 ];
INF_SCATT7                (idx, [1:   4]) = [  2.95502E-04 0.07264  3.24409E-03 0.05176 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  4.01528E-01 0.00035  6.95026E-01 0.00055 ];
INF_SCATTP1               (idx, [1:   4]) = [  1.17236E-01 0.00058  2.04895E-01 0.00154 ];
INF_SCATTP2               (idx, [1:   4]) = [  4.63205E-02 0.00125  7.54482E-02 0.00248 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.34398E-03 0.01569  2.98555E-02 0.00631 ];
INF_SCATTP4               (idx, [1:   4]) = [ -3.66945E-03 0.01417  1.38752E-02 0.01573 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.87399E-04 0.29497  7.38054E-03 0.02039 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.10673E-03 0.01148  4.47094E-03 0.02327 ];
INF_SCATTP7               (idx, [1:   4]) = [  2.95477E-04 0.07274  3.24409E-03 0.05176 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.32696E-01 0.00061  5.02930E-01 0.00101 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.43249E+00 0.00061  6.62785E-01 0.00101 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  3.88449E-03 0.00092  4.94632E-02 0.00132 ];
INF_REMXS                 (idx, [1:   4]) = [  7.63515E-03 0.00167  5.39612E-02 0.00151 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37701E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49395E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.94297E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.92615E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55062E-04 ;
I135_BR                   (idx, 1)        =  8.34918E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  3.97782E-01 0.00036  3.73660E-03 0.00310  4.63241E-03 0.01263  6.90393E-01 0.00050 ];
INF_S1                    (idx, [1:   8]) = [  1.16107E-01 0.00056  1.12901E-03 0.00611  1.61324E-03 0.01803  2.03282E-01 0.00157 ];
INF_S2                    (idx, [1:   8]) = [  4.65577E-02 0.00121 -2.37552E-04 0.02051  8.02154E-04 0.03145  7.46461E-02 0.00273 ];
INF_S3                    (idx, [1:   8]) = [  4.74415E-03 0.01484 -3.99877E-04 0.01024  2.64460E-04 0.08929  2.95910E-02 0.00622 ];
INF_S4                    (idx, [1:   8]) = [ -3.51119E-03 0.01509 -1.58051E-04 0.02006  2.88178E-05 0.69825  1.38464E-02 0.01587 ];
INF_S5                    (idx, [1:   8]) = [  2.15609E-04 0.25237 -2.82703E-05 0.09730 -8.96039E-05 0.25295  7.47014E-03 0.02270 ];
INF_S6                    (idx, [1:   8]) = [  2.11512E-03 0.01133 -8.15598E-06 0.05108 -1.09180E-04 0.11291  4.58012E-03 0.02101 ];
INF_S7                    (idx, [1:   8]) = [  3.07779E-04 0.07196 -1.22771E-05 0.10702 -9.55892E-05 0.14943  3.33968E-03 0.04705 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.97792E-01 0.00036  3.73660E-03 0.00310  4.63241E-03 0.01263  6.90393E-01 0.00050 ];
INF_SP1                   (idx, [1:   8]) = [  1.16107E-01 0.00056  1.12901E-03 0.00611  1.61324E-03 0.01803  2.03282E-01 0.00157 ];
INF_SP2                   (idx, [1:   8]) = [  4.65581E-02 0.00121 -2.37552E-04 0.02051  8.02154E-04 0.03145  7.46461E-02 0.00273 ];
INF_SP3                   (idx, [1:   8]) = [  4.74386E-03 0.01483 -3.99877E-04 0.01024  2.64460E-04 0.08929  2.95910E-02 0.00622 ];
INF_SP4                   (idx, [1:   8]) = [ -3.51140E-03 0.01511 -1.58051E-04 0.02006  2.88178E-05 0.69825  1.38464E-02 0.01587 ];
INF_SP5                   (idx, [1:   8]) = [  2.15669E-04 0.25249 -2.82703E-05 0.09730 -8.96039E-05 0.25295  7.47014E-03 0.02270 ];
INF_SP6                   (idx, [1:   8]) = [  2.11488E-03 0.01132 -8.15598E-06 0.05108 -1.09180E-04 0.11291  4.58012E-03 0.02101 ];
INF_SP7                   (idx, [1:   8]) = [  3.07754E-04 0.07207 -1.22771E-05 0.10702 -9.55892E-05 0.14943  3.33968E-03 0.04705 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  1.00211E-01 0.00121  2.23128E-01 0.00776 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  1.31129E-01 0.00147  4.31950E-01 0.01100 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  1.11849E-01 0.00137  2.56717E-01 0.01356 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  7.47934E-02 0.00181  1.38292E-01 0.01016 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  3.32635E+00 0.00121  1.49427E+00 0.00779 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  2.54204E+00 0.00147  7.72065E-01 0.01092 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  2.98023E+00 0.00137  1.29939E+00 0.01345 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  4.45678E+00 0.00182  2.41135E+00 0.01007 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  9.12183E-03 0.02647  2.69861E-04 0.13532  1.62978E-03 0.05912  1.58369E-03 0.05554  3.41389E-03 0.03812  1.55685E-03 0.05431  6.67767E-04 0.07449 ];
LAMBDA                    (idx, [1:  14]) = [  5.04616E-01 0.02970  1.33594E-02 0.00048  3.25590E-02 0.00067  1.21254E-01 0.00044  3.07737E-01 0.00114  8.67043E-01 0.00137  2.90225E+00 0.00187 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'May 19 2021 21:58:01' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 27])  = 'hexassembly_2D_absorber.inp' ;
WORKING_DIRECTORY         (idx, [1: 26])  = '/nfs/home/jgj/projects/UPV' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode01' ;
CPU_TYPE                  (idx, [1: 46])  = 'AMD Ryzen Threadripper 3970X 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 137367609.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct  4 21:34:37 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct  4 21:35:17 2022' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1664912077112 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 30 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.02283E+00  1.02800E+00  1.02219E+00  1.02498E+00  1.04867E+00  1.05263E+00  1.02384E+00  7.91773E-01  7.91122E-01  1.02798E+00  1.05287E+00  1.02005E+00  7.89558E-01  1.02426E+00  1.05196E+00  1.03010E+00  1.01629E+00  1.04779E+00  1.02655E+00  1.05345E+00  7.90504E-01  1.05409E+00  1.02574E+00  1.01319E+00  1.01636E+00  1.02654E+00  1.02312E+00  1.02609E+00  1.02508E+00  1.05239E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 35])  = '/opt/xsdata/xsdir_files/xsdata_endf' ;
DECAY_DATA_FILE_PATH      (idx, [1: 48])  = '/opt/xsdata/JEFF33/Decay_Data/JEFF33-rdd_all.asc' ;
SFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
NFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 3.7E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.10661E-01 0.00051  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  8.89339E-01 6.4E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.79010E-01 0.00014  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  5.25298E-01 0.00012  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.78928E+00 0.00027  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  1.68602E+01 0.00035  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  1.60739E+01 0.00037  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.45258E+01 0.00051  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.36088E+00 0.00035  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 9997390 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.49683E+01 ;
RUNNING_TIME              (idx, 1)        =  6.65317E-01 ;
INIT_TIME                 (idx, [1:  2])  = [  4.35167E-02  4.35167E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.68333E-03  1.68333E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  6.19433E-01  6.19433E-01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  6.64300E-01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 22.49806 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.41043E+01 0.00057 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  7.00794E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128701.35 ;
ALLOC_MEMSIZE             (idx, 1)        = 1400.32;
MEMSIZE                   (idx, 1)        = 1091.76;
XS_MEMSIZE                (idx, 1)        = 387.35;
MAT_MEMSIZE               (idx, 1)        = 27.72;
RES_MEMSIZE               (idx, 1)        = 8.37;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 668.32;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 308.56;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 23 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 241498 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 2 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 14 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 14 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 309 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.99831E-06 0.00013  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.62158E-01 0.00173 ];
U235_FISS                 (idx, [1:   4]) = [  9.00204E-02 0.00098  9.12290E-01 0.00029 ];
U238_FISS                 (idx, [1:   4]) = [  8.65552E-03 0.00345  8.77097E-02 0.00305 ];
U235_CAPT                 (idx, [1:   4]) = [  2.86726E-02 0.00179  2.49119E-01 0.00155 ];
U238_CAPT                 (idx, [1:   4]) = [  4.29957E-02 0.00150  3.73560E-01 0.00112 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 9997390 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 3.65920E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 1150434 1.15117E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 986249 9.86921E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 7860707 7.86557E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -8.17701E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  3.20185E-12 0.00076 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  2.44351E-01 0.00075 ];
TOT_FISSRATE              (idx, [1:   2]) = [  9.86052E-02 0.00076 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  1.14971E-01 0.00061 ];
TOT_ABSRATE               (idx, [1:   2]) = [  2.13577E-01 0.00060 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.99831E-01 0.00013 ];
TOT_FLUX                  (idx, [1:   2]) = [  3.83988E+01 0.00029 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.86423E-01 0.00016 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  1.60777E+01 0.00034 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.86402E+00 0.00073 ];
SIX_FF_F                  (idx, [1:   2]) = [  8.52868E-01 0.00049 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.82670E-01 0.00110 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.55022E+00 0.00120 ];
SIX_FF_LF                 (idx, [1:   2]) = [  2.80565E-01 0.00059 ];
SIX_FF_LT                 (idx, [1:   2]) = [  7.60762E-01 0.00034 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14585E+00 0.00080 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  2.44572E-01 0.00097 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.47807E+00 1.6E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02671E+02 1.7E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  2.44561E-01 0.00098  2.42443E-01 0.00097  2.12865E-03 0.01050 ];
IMP_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
COL_KEFF                  (idx, [1:   2]) = [  2.44393E-01 0.00080 ];
ABS_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14599E+00 0.00034 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.31676E+01 0.00048 ];
IMP_ALF                   (idx, [1:   2]) = [  1.31807E+01 0.00027 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  3.83063E-05 0.00631 ];
IMP_EALF                  (idx, [1:   2]) = [  3.77563E-05 0.00355 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  3.83173E-01 0.00272 ];
IMP_AFGE                  (idx, [1:   2]) = [  3.82479E-01 0.00091 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  3.01927E-02 0.00419  9.24995E-04 0.02176  5.03528E-03 0.00845  4.99499E-03 0.00923  1.16770E-02 0.00547  5.34205E-03 0.00867  2.21837E-03 0.01298 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  5.11675E-01 0.00488  1.33611E-02 0.00014  3.25403E-02 0.00013  1.21199E-01 7.0E-05  3.07058E-01 0.00018  8.66487E-01 0.00026  2.90811E+00 0.00043 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  8.81649E-03 0.00816  2.72268E-04 0.05149  1.50312E-03 0.01995  1.45803E-03 0.02070  3.35929E-03 0.01468  1.58926E-03 0.01893  6.34525E-04 0.03381 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  5.08559E-01 0.01310  1.33544E-02 0.00022  3.25614E-02 0.00032  1.21174E-01 0.00017  3.06933E-01 0.00041  8.66936E-01 0.00062  2.90992E+00 0.00097 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  4.85774E-05 0.00236  4.85134E-05 0.00235  5.61001E-05 0.02312 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.18788E-05 0.00207  1.18632E-05 0.00207  1.37153E-05 0.02293 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  8.68974E-03 0.01043  2.70582E-04 0.06334  1.47407E-03 0.02434  1.47764E-03 0.02685  3.28101E-03 0.01582  1.55596E-03 0.02523  6.30468E-04 0.04396 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  5.07907E-01 0.01635  1.33534E-02 0.00034  3.25553E-02 0.00040  1.21201E-01 0.00021  3.07118E-01 0.00056  8.67255E-01 0.00084  2.90641E+00 0.00138 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.85141E-05 0.00658  4.84099E-05 0.00652  5.75531E-05 0.05433 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.18632E-05 0.00647  1.18378E-05 0.00641  1.40664E-05 0.05415 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  8.66956E-03 0.04306  4.20749E-04 0.22546  1.41031E-03 0.09430  1.20524E-03 0.09711  3.31507E-03 0.07513  1.38061E-03 0.11194  9.37594E-04 0.14050 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  5.79172E-01 0.06119  1.33570E-02 0.00109  3.25956E-02 0.00116  1.21224E-01 0.00076  3.06904E-01 0.00173  8.66676E-01 0.00283  2.90741E+00 0.00357 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  8.72988E-03 0.04178  4.15360E-04 0.22586  1.38787E-03 0.09259  1.24556E-03 0.09633  3.32720E-03 0.07154  1.40820E-03 0.11097  9.45694E-04 0.13655 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  5.86308E-01 0.06010  1.33570E-02 0.00109  3.25985E-02 0.00115  1.21217E-01 0.00074  3.06821E-01 0.00170  8.66841E-01 0.00284  2.90677E+00 0.00356 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.80255E+02 0.04373 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.84821E-05 0.00148 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.18556E-05 0.00104 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  8.55066E-03 0.00690 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.76396E+02 0.00696 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  8.86765E-08 0.00092 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.72310E-06 0.00091  3.72311E-06 0.00093  3.71805E-06 0.00789 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.86564E-05 0.00093  1.86590E-05 0.00094  1.84121E-05 0.00995 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  1.28519E-01 0.00095  1.31202E-01 0.00095  4.24282E-02 0.00932 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.01810E+01 0.00891 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  1.60739E+01 0.00037  2.57880E+01 0.00064 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '4' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  6.36669E+04 0.00923  3.09356E+05 0.00285  7.02396E+05 0.00366  9.66053E+05 0.00450  9.46924E+05 0.00348  9.41620E+05 0.00277  5.60601E+05 0.00264  4.78390E+05 0.00268  5.56867E+05 0.00262  4.41331E+05 0.00337  3.61984E+05 0.00217  2.94531E+05 0.00256  2.59601E+05 0.00355  2.12510E+05 0.00400  1.95761E+05 0.00211  1.35837E+05 0.00245  6.78930E+04 0.00215  2.09751E+05 0.00126  1.72278E+05 0.00146  2.90165E+05 0.00347  2.57469E+05 0.00126  1.67381E+05 0.00335  9.12663E+04 0.00322  9.46736E+04 0.00493  8.19618E+04 0.00456  6.82487E+04 0.00566  1.02942E+05 0.00194  2.25154E+04 0.00364  2.70092E+04 0.00257  2.43445E+04 0.00401  1.34697E+04 0.00883  2.32191E+04 0.00384  1.52185E+04 0.00916  1.22468E+04 0.00715  2.19946E+03 0.01136  2.10897E+03 0.01622  2.17648E+03 0.02118  2.29530E+03 0.01996  2.25664E+03 0.02204  2.11587E+03 0.01121  2.30521E+03 0.01612  2.10149E+03 0.00857  3.96466E+03 0.02323  6.28013E+03 0.01150  7.92194E+03 0.00720  2.04389E+04 0.00332  2.10363E+04 0.00497  2.22815E+04 0.00582  1.41937E+04 0.00804  9.84358E+03 0.01260  7.25421E+03 0.00718  7.87216E+03 0.00877  1.36331E+04 0.00585  1.59639E+04 0.00496  2.59788E+04 0.00602  3.11150E+04 0.00462  3.54059E+04 0.00449  1.84355E+04 0.00769  1.17592E+04 0.00514  7.81108E+03 0.00307  6.55412E+03 0.00344  6.05198E+03 0.00524  4.85187E+03 0.00227  3.11736E+03 0.01692  2.83074E+03 0.00965  2.42565E+03 0.00728  1.99128E+03 0.01232  1.46037E+03 0.00784  9.16074E+02 0.01081  3.32950E+02 0.01328 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.19311E+00 0.00125 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.61205E+00 0.00256  1.36535E-01 0.00277 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.06447E-01 0.00025  7.40392E-01 0.00117 ];
INF_CAPT                  (idx, [1:   4]) = [  2.24174E-03 0.00125  1.69120E-02 0.00227 ];
INF_ABS                   (idx, [1:   4]) = [  3.82151E-03 0.00103  4.98611E-02 0.00439 ];
INF_FISS                  (idx, [1:   4]) = [  1.57977E-03 0.00083  3.29491E-02 0.00548 ];
INF_NSF                   (idx, [1:   4]) = [  3.95913E-03 0.00081  8.02872E-02 0.00548 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.50614E+00 2.9E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02943E+02 6.7E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  2.65165E-08 0.00125  1.96097E-06 0.00085 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  4.02606E-01 0.00024  6.90431E-01 0.00146 ];
INF_SCATT1                (idx, [1:   4]) = [  1.17705E-01 0.00034  2.03476E-01 0.00251 ];
INF_SCATT2                (idx, [1:   4]) = [  4.63648E-02 0.00115  7.48180E-02 0.00165 ];
INF_SCATT3                (idx, [1:   4]) = [  4.17159E-03 0.00904  3.00223E-02 0.00632 ];
INF_SCATT4                (idx, [1:   4]) = [ -3.78030E-03 0.00664  1.34885E-02 0.01212 ];
INF_SCATT5                (idx, [1:   4]) = [  1.16025E-04 0.13525  7.32319E-03 0.02516 ];
INF_SCATT6                (idx, [1:   4]) = [  2.14880E-03 0.00669  4.72381E-03 0.05128 ];
INF_SCATT7                (idx, [1:   4]) = [  2.73146E-04 0.09747  3.16081E-03 0.08801 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  4.02616E-01 0.00024  6.90431E-01 0.00146 ];
INF_SCATTP1               (idx, [1:   4]) = [  1.17705E-01 0.00034  2.03476E-01 0.00251 ];
INF_SCATTP2               (idx, [1:   4]) = [  4.63646E-02 0.00115  7.48180E-02 0.00165 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.17147E-03 0.00904  3.00223E-02 0.00632 ];
INF_SCATTP4               (idx, [1:   4]) = [ -3.78039E-03 0.00662  1.34885E-02 0.01212 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.16137E-04 0.13520  7.32319E-03 0.02516 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.14875E-03 0.00667  4.72381E-03 0.05128 ];
INF_SCATTP7               (idx, [1:   4]) = [  2.73192E-04 0.09747  3.16081E-03 0.08801 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.32729E-01 0.00025  5.00650E-01 0.00081 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.43228E+00 0.00025  6.65803E-01 0.00081 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  3.81195E-03 0.00107  4.98611E-02 0.00439 ];
INF_REMXS                 (idx, [1:   4]) = [  7.21993E-03 0.00213  5.47746E-02 0.00351 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37701E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49395E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.94297E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.92615E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55062E-04 ;
I135_BR                   (idx, 1)        =  8.34918E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  3.99228E-01 0.00022  3.37893E-03 0.00320  4.81362E-03 0.01200  6.85617E-01 0.00146 ];
INF_S1                    (idx, [1:   8]) = [  1.16694E-01 0.00032  1.01116E-03 0.00480  1.65766E-03 0.01366  2.01818E-01 0.00253 ];
INF_S2                    (idx, [1:   8]) = [  4.65976E-02 0.00109 -2.32778E-04 0.01411  8.12606E-04 0.02503  7.40054E-02 0.00164 ];
INF_S3                    (idx, [1:   8]) = [  4.53721E-03 0.00805 -3.65617E-04 0.00497  2.70198E-04 0.10688  2.97521E-02 0.00687 ];
INF_S4                    (idx, [1:   8]) = [ -3.63638E-03 0.00678 -1.43915E-04 0.00983 -2.30271E-06 1.00000  1.34908E-02 0.01129 ];
INF_S5                    (idx, [1:   8]) = [  1.38056E-04 0.11051 -2.20316E-05 0.07619 -1.01975E-04 0.19595  7.42516E-03 0.02420 ];
INF_S6                    (idx, [1:   8]) = [  2.15456E-03 0.00692 -5.76047E-06 0.32199 -1.04411E-04 0.20882  4.82822E-03 0.04902 ];
INF_S7                    (idx, [1:   8]) = [  2.79975E-04 0.09057 -6.82914E-06 0.37557 -8.72763E-05 0.18081  3.24808E-03 0.08185 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.99237E-01 0.00022  3.37893E-03 0.00320  4.81362E-03 0.01200  6.85617E-01 0.00146 ];
INF_SP1                   (idx, [1:   8]) = [  1.16694E-01 0.00032  1.01116E-03 0.00480  1.65766E-03 0.01366  2.01818E-01 0.00253 ];
INF_SP2                   (idx, [1:   8]) = [  4.65973E-02 0.00109 -2.32778E-04 0.01411  8.12606E-04 0.02503  7.40054E-02 0.00164 ];
INF_SP3                   (idx, [1:   8]) = [  4.53709E-03 0.00805 -3.65617E-04 0.00497  2.70198E-04 0.10688  2.97521E-02 0.00687 ];
INF_SP4                   (idx, [1:   8]) = [ -3.63648E-03 0.00676 -1.43915E-04 0.00983 -2.30271E-06 1.00000  1.34908E-02 0.01129 ];
INF_SP5                   (idx, [1:   8]) = [  1.38168E-04 0.11050 -2.20316E-05 0.07619 -1.01975E-04 0.19595  7.42516E-03 0.02420 ];
INF_SP6                   (idx, [1:   8]) = [  2.15451E-03 0.00690 -5.76047E-06 0.32199 -1.04411E-04 0.20882  4.82822E-03 0.04902 ];
INF_SP7                   (idx, [1:   8]) = [  2.80021E-04 0.09056 -6.82914E-06 0.37557 -8.72763E-05 0.18081  3.24808E-03 0.08185 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  7.46950E-02 0.00266  1.42972E-01 0.00592 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  8.83628E-02 0.00252  2.55158E-01 0.00915 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  8.74266E-02 0.00371  1.62926E-01 0.01044 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  5.74448E-02 0.00229  9.15713E-02 0.01045 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  4.46272E+00 0.00265  2.33179E+00 0.00597 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  3.77242E+00 0.00252  1.30682E+00 0.00912 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  3.81293E+00 0.00372  2.04681E+00 0.01045 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  5.80280E+00 0.00229  3.64174E+00 0.01047 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  8.54620E-03 0.02701  2.31257E-04 0.13208  1.32918E-03 0.05141  1.48996E-03 0.05340  3.18756E-03 0.04924  1.67554E-03 0.07109  6.32701E-04 0.08671 ];
LAMBDA                    (idx, [1:  14]) = [  5.22828E-01 0.03664  1.33526E-02 0.00042  3.25516E-02 0.00077  1.21156E-01 0.00038  3.07120E-01 0.00094  8.67299E-01 0.00153  2.91029E+00 0.00196 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'May 19 2021 21:58:01' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 27])  = 'hexassembly_2D_absorber.inp' ;
WORKING_DIRECTORY         (idx, [1: 26])  = '/nfs/home/jgj/projects/UPV' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode01' ;
CPU_TYPE                  (idx, [1: 46])  = 'AMD Ryzen Threadripper 3970X 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 137367609.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct  4 21:34:37 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct  4 21:35:17 2022' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1664912077112 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 30 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.02283E+00  1.02800E+00  1.02219E+00  1.02498E+00  1.04867E+00  1.05263E+00  1.02384E+00  7.91773E-01  7.91122E-01  1.02798E+00  1.05287E+00  1.02005E+00  7.89558E-01  1.02426E+00  1.05196E+00  1.03010E+00  1.01629E+00  1.04779E+00  1.02655E+00  1.05345E+00  7.90504E-01  1.05409E+00  1.02574E+00  1.01319E+00  1.01636E+00  1.02654E+00  1.02312E+00  1.02609E+00  1.02508E+00  1.05239E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 35])  = '/opt/xsdata/xsdir_files/xsdata_endf' ;
DECAY_DATA_FILE_PATH      (idx, [1: 48])  = '/opt/xsdata/JEFF33/Decay_Data/JEFF33-rdd_all.asc' ;
SFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
NFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 3.7E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.10661E-01 0.00051  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  8.89339E-01 6.4E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.79010E-01 0.00014  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  5.25298E-01 0.00012  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.78928E+00 0.00027  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  1.68602E+01 0.00035  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  1.60739E+01 0.00037  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.45258E+01 0.00051  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.36088E+00 0.00035  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 9997390 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.49683E+01 ;
RUNNING_TIME              (idx, 1)        =  6.65317E-01 ;
INIT_TIME                 (idx, [1:  2])  = [  4.35167E-02  4.35167E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.68333E-03  1.68333E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  6.19433E-01  6.19433E-01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  6.64300E-01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 22.49806 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.41043E+01 0.00057 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  7.00794E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128701.35 ;
ALLOC_MEMSIZE             (idx, 1)        = 1400.32;
MEMSIZE                   (idx, 1)        = 1091.76;
XS_MEMSIZE                (idx, 1)        = 387.35;
MAT_MEMSIZE               (idx, 1)        = 27.72;
RES_MEMSIZE               (idx, 1)        = 8.37;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 668.32;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 308.56;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 23 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 241498 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 2 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 14 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 14 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 309 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.99831E-06 0.00013  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.62158E-01 0.00173 ];
U235_FISS                 (idx, [1:   4]) = [  9.00204E-02 0.00098  9.12290E-01 0.00029 ];
U238_FISS                 (idx, [1:   4]) = [  8.65552E-03 0.00345  8.77097E-02 0.00305 ];
U235_CAPT                 (idx, [1:   4]) = [  2.86726E-02 0.00179  2.49119E-01 0.00155 ];
U238_CAPT                 (idx, [1:   4]) = [  4.29957E-02 0.00150  3.73560E-01 0.00112 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 9997390 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 3.65920E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 1150434 1.15117E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 986249 9.86921E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 7860707 7.86557E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -8.17701E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  3.20185E-12 0.00076 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  2.44351E-01 0.00075 ];
TOT_FISSRATE              (idx, [1:   2]) = [  9.86052E-02 0.00076 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  1.14971E-01 0.00061 ];
TOT_ABSRATE               (idx, [1:   2]) = [  2.13577E-01 0.00060 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.99831E-01 0.00013 ];
TOT_FLUX                  (idx, [1:   2]) = [  3.83988E+01 0.00029 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.86423E-01 0.00016 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  1.60777E+01 0.00034 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.86402E+00 0.00073 ];
SIX_FF_F                  (idx, [1:   2]) = [  8.52868E-01 0.00049 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.82670E-01 0.00110 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.55022E+00 0.00120 ];
SIX_FF_LF                 (idx, [1:   2]) = [  2.80565E-01 0.00059 ];
SIX_FF_LT                 (idx, [1:   2]) = [  7.60762E-01 0.00034 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14585E+00 0.00080 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  2.44572E-01 0.00097 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.47807E+00 1.6E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02671E+02 1.7E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  2.44561E-01 0.00098  2.42443E-01 0.00097  2.12865E-03 0.01050 ];
IMP_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
COL_KEFF                  (idx, [1:   2]) = [  2.44393E-01 0.00080 ];
ABS_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14599E+00 0.00034 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.31676E+01 0.00048 ];
IMP_ALF                   (idx, [1:   2]) = [  1.31807E+01 0.00027 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  3.83063E-05 0.00631 ];
IMP_EALF                  (idx, [1:   2]) = [  3.77563E-05 0.00355 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  3.83173E-01 0.00272 ];
IMP_AFGE                  (idx, [1:   2]) = [  3.82479E-01 0.00091 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  3.01927E-02 0.00419  9.24995E-04 0.02176  5.03528E-03 0.00845  4.99499E-03 0.00923  1.16770E-02 0.00547  5.34205E-03 0.00867  2.21837E-03 0.01298 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  5.11675E-01 0.00488  1.33611E-02 0.00014  3.25403E-02 0.00013  1.21199E-01 7.0E-05  3.07058E-01 0.00018  8.66487E-01 0.00026  2.90811E+00 0.00043 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  8.81649E-03 0.00816  2.72268E-04 0.05149  1.50312E-03 0.01995  1.45803E-03 0.02070  3.35929E-03 0.01468  1.58926E-03 0.01893  6.34525E-04 0.03381 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  5.08559E-01 0.01310  1.33544E-02 0.00022  3.25614E-02 0.00032  1.21174E-01 0.00017  3.06933E-01 0.00041  8.66936E-01 0.00062  2.90992E+00 0.00097 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  4.85774E-05 0.00236  4.85134E-05 0.00235  5.61001E-05 0.02312 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.18788E-05 0.00207  1.18632E-05 0.00207  1.37153E-05 0.02293 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  8.68974E-03 0.01043  2.70582E-04 0.06334  1.47407E-03 0.02434  1.47764E-03 0.02685  3.28101E-03 0.01582  1.55596E-03 0.02523  6.30468E-04 0.04396 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  5.07907E-01 0.01635  1.33534E-02 0.00034  3.25553E-02 0.00040  1.21201E-01 0.00021  3.07118E-01 0.00056  8.67255E-01 0.00084  2.90641E+00 0.00138 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.85141E-05 0.00658  4.84099E-05 0.00652  5.75531E-05 0.05433 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.18632E-05 0.00647  1.18378E-05 0.00641  1.40664E-05 0.05415 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  8.66956E-03 0.04306  4.20749E-04 0.22546  1.41031E-03 0.09430  1.20524E-03 0.09711  3.31507E-03 0.07513  1.38061E-03 0.11194  9.37594E-04 0.14050 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  5.79172E-01 0.06119  1.33570E-02 0.00109  3.25956E-02 0.00116  1.21224E-01 0.00076  3.06904E-01 0.00173  8.66676E-01 0.00283  2.90741E+00 0.00357 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  8.72988E-03 0.04178  4.15360E-04 0.22586  1.38787E-03 0.09259  1.24556E-03 0.09633  3.32720E-03 0.07154  1.40820E-03 0.11097  9.45694E-04 0.13655 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  5.86308E-01 0.06010  1.33570E-02 0.00109  3.25985E-02 0.00115  1.21217E-01 0.00074  3.06821E-01 0.00170  8.66841E-01 0.00284  2.90677E+00 0.00356 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.80255E+02 0.04373 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.84821E-05 0.00148 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.18556E-05 0.00104 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  8.55066E-03 0.00690 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.76396E+02 0.00696 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  8.86765E-08 0.00092 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.72310E-06 0.00091  3.72311E-06 0.00093  3.71805E-06 0.00789 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.86564E-05 0.00093  1.86590E-05 0.00094  1.84121E-05 0.00995 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  1.28519E-01 0.00095  1.31202E-01 0.00095  4.24282E-02 0.00932 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.01810E+01 0.00891 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  1.60739E+01 0.00037  2.57880E+01 0.00064 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '5' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  6.38217E+04 0.00590  3.06658E+05 0.00593  7.00391E+05 0.00528  9.64053E+05 0.00374  9.45536E+05 0.00443  9.38904E+05 0.00349  5.58696E+05 0.00321  4.75832E+05 0.00370  5.56676E+05 0.00302  4.40275E+05 0.00335  3.59965E+05 0.00269  2.93781E+05 0.00204  2.58483E+05 0.00328  2.12048E+05 0.00500  1.94388E+05 0.00460  1.35904E+05 0.00398  6.79559E+04 0.00215  2.09306E+05 0.00447  1.72578E+05 0.00353  2.91888E+05 0.00339  2.57430E+05 0.00439  1.68358E+05 0.00476  9.12123E+04 0.00239  9.53485E+04 0.00408  8.16619E+04 0.00452  6.74793E+04 0.00415  1.02556E+05 0.00482  2.21085E+04 0.00687  2.70419E+04 0.01009  2.42759E+04 0.00481  1.34592E+04 0.01045  2.30037E+04 0.00969  1.50147E+04 0.01070  1.20818E+04 0.00646  2.19573E+03 0.00558  2.15437E+03 0.01435  2.18960E+03 0.01715  2.21394E+03 0.02553  2.21351E+03 0.01898  2.19285E+03 0.01544  2.24480E+03 0.00845  2.12309E+03 0.01463  3.96187E+03 0.00663  6.27832E+03 0.01238  7.90491E+03 0.00700  2.04400E+04 0.00345  2.10780E+04 0.00499  2.19415E+04 0.00369  1.39488E+04 0.00917  9.74605E+03 0.00638  7.19423E+03 0.01183  7.86513E+03 0.00312  1.34636E+04 0.00468  1.59389E+04 0.00488  2.58668E+04 0.00465  3.09502E+04 0.00754  3.51862E+04 0.00316  1.82367E+04 0.00523  1.16780E+04 0.00197  7.73309E+03 0.00994  6.48294E+03 0.00465  6.01575E+03 0.00333  4.81396E+03 0.00827  3.11183E+03 0.00498  2.80104E+03 0.01280  2.42749E+03 0.00690  2.00236E+03 0.00390  1.47772E+03 0.01048  9.52045E+02 0.01031  3.28732E+02 0.00755 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.19448E+00 0.00169 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.60136E+00 0.00313  1.35598E-01 0.00278 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.06645E-01 0.00043  7.41261E-01 0.00037 ];
INF_CAPT                  (idx, [1:   4]) = [  2.25414E-03 0.00178  1.68825E-02 0.00082 ];
INF_ABS                   (idx, [1:   4]) = [  3.83974E-03 0.00144  4.96512E-02 0.00151 ];
INF_FISS                  (idx, [1:   4]) = [  1.58560E-03 0.00123  3.27687E-02 0.00187 ];
INF_NSF                   (idx, [1:   4]) = [  3.97313E-03 0.00125  7.98476E-02 0.00187 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.50575E+00 6.0E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02939E+02 7.0E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  2.64966E-08 0.00218  1.96285E-06 0.00046 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  4.02802E-01 0.00045  6.91629E-01 0.00052 ];
INF_SCATT1                (idx, [1:   4]) = [  1.17933E-01 0.00088  2.03372E-01 0.00247 ];
INF_SCATT2                (idx, [1:   4]) = [  4.65140E-02 0.00117  7.53173E-02 0.00339 ];
INF_SCATT3                (idx, [1:   4]) = [  4.21930E-03 0.00708  3.00120E-02 0.00684 ];
INF_SCATT4                (idx, [1:   4]) = [ -3.75754E-03 0.00330  1.38950E-02 0.01700 ];
INF_SCATT5                (idx, [1:   4]) = [  1.50416E-04 0.28557  7.45969E-03 0.02001 ];
INF_SCATT6                (idx, [1:   4]) = [  2.12371E-03 0.01174  4.60217E-03 0.02919 ];
INF_SCATT7                (idx, [1:   4]) = [  2.67226E-04 0.07954  3.02405E-03 0.07777 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  4.02812E-01 0.00045  6.91629E-01 0.00052 ];
INF_SCATTP1               (idx, [1:   4]) = [  1.17933E-01 0.00088  2.03372E-01 0.00247 ];
INF_SCATTP2               (idx, [1:   4]) = [  4.65141E-02 0.00117  7.53173E-02 0.00339 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.21932E-03 0.00705  3.00120E-02 0.00684 ];
INF_SCATTP4               (idx, [1:   4]) = [ -3.75736E-03 0.00331  1.38950E-02 0.01700 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.50366E-04 0.28521  7.45969E-03 0.02001 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.12381E-03 0.01168  4.60217E-03 0.02919 ];
INF_SCATTP7               (idx, [1:   4]) = [  2.67248E-04 0.07952  3.02405E-03 0.07777 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.32553E-01 0.00036  5.01386E-01 0.00091 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.43337E+00 0.00036  6.64826E-01 0.00091 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  3.82994E-03 0.00150  4.96512E-02 0.00151 ];
INF_REMXS                 (idx, [1:   4]) = [  7.22966E-03 0.00224  5.45852E-02 0.00396 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37701E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49395E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.94297E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.92615E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55062E-04 ;
I135_BR                   (idx, 1)        =  8.34918E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  3.99415E-01 0.00044  3.38624E-03 0.00375  4.95300E-03 0.01478  6.86676E-01 0.00060 ];
INF_S1                    (idx, [1:   8]) = [  1.16917E-01 0.00088  1.01595E-03 0.00543  1.69613E-03 0.04253  2.01676E-01 0.00268 ];
INF_S2                    (idx, [1:   8]) = [  4.67489E-02 0.00112 -2.34930E-04 0.02023  8.54770E-04 0.02815  7.44626E-02 0.00313 ];
INF_S3                    (idx, [1:   8]) = [  4.58259E-03 0.00664 -3.63292E-04 0.00844  2.94023E-04 0.07825  2.97179E-02 0.00689 ];
INF_S4                    (idx, [1:   8]) = [ -3.62069E-03 0.00320 -1.36847E-04 0.01071  4.89482E-05 0.27435  1.38461E-02 0.01710 ];
INF_S5                    (idx, [1:   8]) = [  1.67546E-04 0.25269 -1.71294E-05 0.15848 -4.48077E-05 0.15681  7.50450E-03 0.02026 ];
INF_S6                    (idx, [1:   8]) = [  2.13010E-03 0.01191 -6.38882E-06 0.31539 -7.23059E-05 0.24762  4.67448E-03 0.02927 ];
INF_S7                    (idx, [1:   8]) = [  2.76087E-04 0.07844 -8.86073E-06 0.14715 -7.34179E-05 0.09580  3.09746E-03 0.07614 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.99425E-01 0.00044  3.38624E-03 0.00375  4.95300E-03 0.01478  6.86676E-01 0.00060 ];
INF_SP1                   (idx, [1:   8]) = [  1.16917E-01 0.00088  1.01595E-03 0.00543  1.69613E-03 0.04253  2.01676E-01 0.00268 ];
INF_SP2                   (idx, [1:   8]) = [  4.67490E-02 0.00112 -2.34930E-04 0.02023  8.54770E-04 0.02815  7.44626E-02 0.00313 ];
INF_SP3                   (idx, [1:   8]) = [  4.58261E-03 0.00661 -3.63292E-04 0.00844  2.94023E-04 0.07825  2.97179E-02 0.00689 ];
INF_SP4                   (idx, [1:   8]) = [ -3.62052E-03 0.00320 -1.36847E-04 0.01071  4.89482E-05 0.27435  1.38461E-02 0.01710 ];
INF_SP5                   (idx, [1:   8]) = [  1.67496E-04 0.25239 -1.71294E-05 0.15848 -4.48077E-05 0.15681  7.50450E-03 0.02026 ];
INF_SP6                   (idx, [1:   8]) = [  2.13020E-03 0.01185 -6.38882E-06 0.31539 -7.23059E-05 0.24762  4.67448E-03 0.02927 ];
INF_SP7                   (idx, [1:   8]) = [  2.76109E-04 0.07841 -8.86073E-06 0.14715 -7.34179E-05 0.09580  3.09746E-03 0.07614 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  6.07303E-02 0.00308  1.17870E-01 0.00415 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  6.77391E-02 0.00393  1.98825E-01 0.01010 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  7.32852E-02 0.00291  1.33570E-01 0.00976 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  4.76400E-02 0.00274  7.73319E-02 0.00471 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  5.48896E+00 0.00308  2.82817E+00 0.00418 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  4.92115E+00 0.00393  1.67720E+00 0.01005 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  4.54859E+00 0.00292  2.49651E+00 0.00963 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  6.99712E+00 0.00273  4.31081E+00 0.00471 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  8.76114E-03 0.02309  3.32942E-04 0.14635  1.45270E-03 0.05727  1.50814E-03 0.05461  3.21038E-03 0.03886  1.55565E-03 0.06024  7.01327E-04 0.08937 ];
LAMBDA                    (idx, [1:  14]) = [  5.23840E-01 0.03628  1.33598E-02 0.00045  3.25485E-02 0.00080  1.21209E-01 0.00037  3.06573E-01 0.00087  8.67030E-01 0.00153  2.90778E+00 0.00189 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'May 19 2021 21:58:01' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 27])  = 'hexassembly_2D_absorber.inp' ;
WORKING_DIRECTORY         (idx, [1: 26])  = '/nfs/home/jgj/projects/UPV' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode01' ;
CPU_TYPE                  (idx, [1: 46])  = 'AMD Ryzen Threadripper 3970X 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 137367609.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct  4 21:34:37 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct  4 21:35:17 2022' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1664912077112 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 30 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.02283E+00  1.02800E+00  1.02219E+00  1.02498E+00  1.04867E+00  1.05263E+00  1.02384E+00  7.91773E-01  7.91122E-01  1.02798E+00  1.05287E+00  1.02005E+00  7.89558E-01  1.02426E+00  1.05196E+00  1.03010E+00  1.01629E+00  1.04779E+00  1.02655E+00  1.05345E+00  7.90504E-01  1.05409E+00  1.02574E+00  1.01319E+00  1.01636E+00  1.02654E+00  1.02312E+00  1.02609E+00  1.02508E+00  1.05239E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 35])  = '/opt/xsdata/xsdir_files/xsdata_endf' ;
DECAY_DATA_FILE_PATH      (idx, [1: 48])  = '/opt/xsdata/JEFF33/Decay_Data/JEFF33-rdd_all.asc' ;
SFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
NFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 3.7E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.10661E-01 0.00051  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  8.89339E-01 6.4E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.79010E-01 0.00014  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  5.25298E-01 0.00012  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.78928E+00 0.00027  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  1.68602E+01 0.00035  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  1.60739E+01 0.00037  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.45258E+01 0.00051  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.36088E+00 0.00035  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 9997390 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.49683E+01 ;
RUNNING_TIME              (idx, 1)        =  6.65333E-01 ;
INIT_TIME                 (idx, [1:  2])  = [  4.35167E-02  4.35167E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.68333E-03  1.68333E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  6.19433E-01  6.19433E-01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  6.64300E-01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 22.49750 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.41043E+01 0.00057 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  7.00777E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128701.35 ;
ALLOC_MEMSIZE             (idx, 1)        = 1400.32;
MEMSIZE                   (idx, 1)        = 1091.76;
XS_MEMSIZE                (idx, 1)        = 387.35;
MAT_MEMSIZE               (idx, 1)        = 27.72;
RES_MEMSIZE               (idx, 1)        = 8.37;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 668.32;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 308.56;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 23 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 241498 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 2 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 14 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 14 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 309 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.99831E-06 0.00013  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.62158E-01 0.00173 ];
U235_FISS                 (idx, [1:   4]) = [  9.00204E-02 0.00098  9.12290E-01 0.00029 ];
U238_FISS                 (idx, [1:   4]) = [  8.65552E-03 0.00345  8.77097E-02 0.00305 ];
U235_CAPT                 (idx, [1:   4]) = [  2.86726E-02 0.00179  2.49119E-01 0.00155 ];
U238_CAPT                 (idx, [1:   4]) = [  4.29957E-02 0.00150  3.73560E-01 0.00112 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 9997390 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 3.65920E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 1150434 1.15117E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 986249 9.86921E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 7860707 7.86557E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -8.17701E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  3.20185E-12 0.00076 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  2.44351E-01 0.00075 ];
TOT_FISSRATE              (idx, [1:   2]) = [  9.86052E-02 0.00076 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  1.14971E-01 0.00061 ];
TOT_ABSRATE               (idx, [1:   2]) = [  2.13577E-01 0.00060 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.99831E-01 0.00013 ];
TOT_FLUX                  (idx, [1:   2]) = [  3.83988E+01 0.00029 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.86423E-01 0.00016 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  1.60777E+01 0.00034 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.86402E+00 0.00073 ];
SIX_FF_F                  (idx, [1:   2]) = [  8.52868E-01 0.00049 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.82670E-01 0.00110 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.55022E+00 0.00120 ];
SIX_FF_LF                 (idx, [1:   2]) = [  2.80565E-01 0.00059 ];
SIX_FF_LT                 (idx, [1:   2]) = [  7.60762E-01 0.00034 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14585E+00 0.00080 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  2.44572E-01 0.00097 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.47807E+00 1.6E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02671E+02 1.7E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  2.44561E-01 0.00098  2.42443E-01 0.00097  2.12865E-03 0.01050 ];
IMP_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
COL_KEFF                  (idx, [1:   2]) = [  2.44393E-01 0.00080 ];
ABS_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14599E+00 0.00034 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.31676E+01 0.00048 ];
IMP_ALF                   (idx, [1:   2]) = [  1.31807E+01 0.00027 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  3.83063E-05 0.00631 ];
IMP_EALF                  (idx, [1:   2]) = [  3.77563E-05 0.00355 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  3.83173E-01 0.00272 ];
IMP_AFGE                  (idx, [1:   2]) = [  3.82479E-01 0.00091 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  3.01927E-02 0.00419  9.24995E-04 0.02176  5.03528E-03 0.00845  4.99499E-03 0.00923  1.16770E-02 0.00547  5.34205E-03 0.00867  2.21837E-03 0.01298 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  5.11675E-01 0.00488  1.33611E-02 0.00014  3.25403E-02 0.00013  1.21199E-01 7.0E-05  3.07058E-01 0.00018  8.66487E-01 0.00026  2.90811E+00 0.00043 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  8.81649E-03 0.00816  2.72268E-04 0.05149  1.50312E-03 0.01995  1.45803E-03 0.02070  3.35929E-03 0.01468  1.58926E-03 0.01893  6.34525E-04 0.03381 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  5.08559E-01 0.01310  1.33544E-02 0.00022  3.25614E-02 0.00032  1.21174E-01 0.00017  3.06933E-01 0.00041  8.66936E-01 0.00062  2.90992E+00 0.00097 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  4.85774E-05 0.00236  4.85134E-05 0.00235  5.61001E-05 0.02312 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.18788E-05 0.00207  1.18632E-05 0.00207  1.37153E-05 0.02293 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  8.68974E-03 0.01043  2.70582E-04 0.06334  1.47407E-03 0.02434  1.47764E-03 0.02685  3.28101E-03 0.01582  1.55596E-03 0.02523  6.30468E-04 0.04396 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  5.07907E-01 0.01635  1.33534E-02 0.00034  3.25553E-02 0.00040  1.21201E-01 0.00021  3.07118E-01 0.00056  8.67255E-01 0.00084  2.90641E+00 0.00138 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.85141E-05 0.00658  4.84099E-05 0.00652  5.75531E-05 0.05433 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.18632E-05 0.00647  1.18378E-05 0.00641  1.40664E-05 0.05415 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  8.66956E-03 0.04306  4.20749E-04 0.22546  1.41031E-03 0.09430  1.20524E-03 0.09711  3.31507E-03 0.07513  1.38061E-03 0.11194  9.37594E-04 0.14050 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  5.79172E-01 0.06119  1.33570E-02 0.00109  3.25956E-02 0.00116  1.21224E-01 0.00076  3.06904E-01 0.00173  8.66676E-01 0.00283  2.90741E+00 0.00357 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  8.72988E-03 0.04178  4.15360E-04 0.22586  1.38787E-03 0.09259  1.24556E-03 0.09633  3.32720E-03 0.07154  1.40820E-03 0.11097  9.45694E-04 0.13655 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  5.86308E-01 0.06010  1.33570E-02 0.00109  3.25985E-02 0.00115  1.21217E-01 0.00074  3.06821E-01 0.00170  8.66841E-01 0.00284  2.90677E+00 0.00356 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.80255E+02 0.04373 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.84821E-05 0.00148 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.18556E-05 0.00104 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  8.55066E-03 0.00690 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.76396E+02 0.00696 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  8.86765E-08 0.00092 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.72310E-06 0.00091  3.72311E-06 0.00093  3.71805E-06 0.00789 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.86564E-05 0.00093  1.86590E-05 0.00094  1.84121E-05 0.00995 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  1.28519E-01 0.00095  1.31202E-01 0.00095  4.24282E-02 0.00932 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.01810E+01 0.00891 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  1.60739E+01 0.00037  2.57880E+01 0.00064 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '6' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  6.79494E+04 0.00706  3.25092E+05 0.00387  7.41391E+05 0.00248  1.01563E+06 0.00279  9.96880E+05 0.00408  9.84426E+05 0.00241  5.84172E+05 0.00166  4.97642E+05 0.00281  5.79271E+05 0.00203  4.55708E+05 0.00382  3.73224E+05 0.00360  3.02483E+05 0.00410  2.65363E+05 0.00284  2.16338E+05 0.00278  1.97409E+05 0.00371  1.37623E+05 0.00449  6.89071E+04 0.00266  2.12434E+05 0.00372  1.74077E+05 0.00387  2.94481E+05 0.00177  2.60277E+05 0.00388  1.71743E+05 0.00411  9.61130E+04 0.00326  1.01580E+05 0.00573  8.84908E+04 0.00445  7.41278E+04 0.00421  1.14510E+05 0.00417  2.50710E+04 0.00426  3.06026E+04 0.00562  2.76048E+04 0.00664  1.56372E+04 0.00683  2.65460E+04 0.00191  1.72616E+04 0.01087  1.38556E+04 0.00793  2.51722E+03 0.01533  2.49980E+03 0.01039  2.46927E+03 0.01771  2.59994E+03 0.01927  2.58360E+03 0.00769  2.55441E+03 0.00807  2.61599E+03 0.01286  2.46770E+03 0.00646  4.50634E+03 0.01048  7.25382E+03 0.00532  9.13238E+03 0.01046  2.40042E+04 0.00282  2.51290E+04 0.00520  2.68622E+04 0.00375  1.72237E+04 0.00455  1.20956E+04 0.00490  8.98976E+03 0.00989  9.98554E+03 0.00716  1.69153E+04 0.00786  2.00764E+04 0.00786  3.26773E+04 0.00432  3.95173E+04 0.00733  4.53046E+04 0.00580  2.33224E+04 0.00778  1.47818E+04 0.00534  9.85776E+03 0.00666  8.33912E+03 0.00360  7.69891E+03 0.00608  6.15630E+03 0.00310  3.92111E+03 0.01010  3.52751E+03 0.00713  3.06813E+03 0.01273  2.49723E+03 0.00752  1.87707E+03 0.00906  1.19174E+03 0.00700  4.27673E+02 0.01723 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.21655E+00 0.00082 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.80875E+00 0.00209  1.70693E-01 0.00470 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.05103E-01 0.00029  7.43756E-01 0.00052 ];
INF_CAPT                  (idx, [1:   4]) = [  2.28086E-03 0.00148  1.69155E-02 0.00118 ];
INF_ABS                   (idx, [1:   4]) = [  3.88984E-03 0.00098  4.96096E-02 0.00206 ];
INF_FISS                  (idx, [1:   4]) = [  1.60898E-03 0.00087  3.26941E-02 0.00251 ];
INF_NSF                   (idx, [1:   4]) = [  4.03390E-03 0.00087  7.96657E-02 0.00251 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.50711E+00 4.1E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02950E+02 5.0E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  2.79751E-08 0.00059  1.97333E-06 0.00015 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  4.01203E-01 0.00030  6.94433E-01 0.00070 ];
INF_SCATT1                (idx, [1:   4]) = [  1.17079E-01 0.00021  2.04461E-01 0.00210 ];
INF_SCATT2                (idx, [1:   4]) = [  4.62817E-02 0.00069  7.51961E-02 0.00313 ];
INF_SCATT3                (idx, [1:   4]) = [  4.27794E-03 0.00682  2.95297E-02 0.00382 ];
INF_SCATT4                (idx, [1:   4]) = [ -3.67312E-03 0.00665  1.33819E-02 0.01941 ];
INF_SCATT5                (idx, [1:   4]) = [  1.70452E-04 0.20106  7.01365E-03 0.01662 ];
INF_SCATT6                (idx, [1:   4]) = [  2.14371E-03 0.01306  4.54869E-03 0.04059 ];
INF_SCATT7                (idx, [1:   4]) = [  3.10622E-04 0.03941  3.36888E-03 0.01680 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  4.01214E-01 0.00030  6.94433E-01 0.00070 ];
INF_SCATTP1               (idx, [1:   4]) = [  1.17079E-01 0.00021  2.04461E-01 0.00210 ];
INF_SCATTP2               (idx, [1:   4]) = [  4.62820E-02 0.00069  7.51961E-02 0.00313 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.27803E-03 0.00683  2.95297E-02 0.00382 ];
INF_SCATTP4               (idx, [1:   4]) = [ -3.67335E-03 0.00665  1.33819E-02 0.01941 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.70463E-04 0.20099  7.01365E-03 0.01662 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.14388E-03 0.01309  4.54869E-03 0.04059 ];
INF_SCATTP7               (idx, [1:   4]) = [  3.10752E-04 0.03931  3.36888E-03 0.01680 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.32576E-01 0.00040  5.02981E-01 0.00034 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.43323E+00 0.00040  6.62716E-01 0.00034 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  3.87907E-03 0.00108  4.96096E-02 0.00206 ];
INF_REMXS                 (idx, [1:   4]) = [  7.62573E-03 0.00127  5.39089E-02 0.00289 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37701E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49395E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.94297E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.92615E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55062E-04 ;
I135_BR                   (idx, 1)        =  8.34918E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  3.97478E-01 0.00029  3.72536E-03 0.00253  4.58644E-03 0.00642  6.89847E-01 0.00071 ];
INF_S1                    (idx, [1:   8]) = [  1.15958E-01 0.00022  1.12109E-03 0.00260  1.54219E-03 0.01717  2.02919E-01 0.00217 ];
INF_S2                    (idx, [1:   8]) = [  4.65336E-02 0.00068 -2.51901E-04 0.00537  7.22581E-04 0.01517  7.44735E-02 0.00308 ];
INF_S3                    (idx, [1:   8]) = [  4.67951E-03 0.00617 -4.01577E-04 0.00807  2.20940E-04 0.07631  2.93087E-02 0.00408 ];
INF_S4                    (idx, [1:   8]) = [ -3.51346E-03 0.00675 -1.59662E-04 0.02146 -7.02241E-06 1.00000  1.33889E-02 0.01960 ];
INF_S5                    (idx, [1:   8]) = [  1.97949E-04 0.16653 -2.74975E-05 0.08658 -7.62340E-05 0.19228  7.08988E-03 0.01723 ];
INF_S6                    (idx, [1:   8]) = [  2.15185E-03 0.01355 -8.14182E-06 0.27858 -6.98260E-05 0.24041  4.61852E-03 0.03918 ];
INF_S7                    (idx, [1:   8]) = [  3.18865E-04 0.03886 -8.24276E-06 0.20024 -8.29386E-05 0.16501  3.45182E-03 0.01376 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.97488E-01 0.00029  3.72536E-03 0.00253  4.58644E-03 0.00642  6.89847E-01 0.00071 ];
INF_SP1                   (idx, [1:   8]) = [  1.15958E-01 0.00022  1.12109E-03 0.00260  1.54219E-03 0.01717  2.02919E-01 0.00217 ];
INF_SP2                   (idx, [1:   8]) = [  4.65339E-02 0.00068 -2.51901E-04 0.00537  7.22581E-04 0.01517  7.44735E-02 0.00308 ];
INF_SP3                   (idx, [1:   8]) = [  4.67961E-03 0.00618 -4.01577E-04 0.00807  2.20940E-04 0.07631  2.93087E-02 0.00408 ];
INF_SP4                   (idx, [1:   8]) = [ -3.51369E-03 0.00675 -1.59662E-04 0.02146 -7.02241E-06 1.00000  1.33889E-02 0.01960 ];
INF_SP5                   (idx, [1:   8]) = [  1.97961E-04 0.16648 -2.74975E-05 0.08658 -7.62340E-05 0.19228  7.08988E-03 0.01723 ];
INF_SP6                   (idx, [1:   8]) = [  2.15203E-03 0.01358 -8.14182E-06 0.27858 -6.98260E-05 0.24041  4.61852E-03 0.03918 ];
INF_SP7                   (idx, [1:   8]) = [  3.18995E-04 0.03876 -8.24276E-06 0.20024 -8.29386E-05 0.16501  3.45182E-03 0.01376 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  5.34013E-02 0.00211  1.13294E-01 0.00574 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  6.41585E-02 0.00229  2.24920E-01 0.01107 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  5.87103E-02 0.00220  1.01973E-01 0.00900 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  4.24464E-02 0.00216  8.17996E-02 0.00295 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  6.24216E+00 0.00212  2.94258E+00 0.00574 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  5.19557E+00 0.00230  1.48273E+00 0.01097 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  5.67771E+00 0.00221  3.26988E+00 0.00894 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  7.85320E+00 0.00217  4.07514E+00 0.00295 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  8.89127E-03 0.02348  2.35710E-04 0.11258  1.62311E-03 0.06433  1.37079E-03 0.05519  3.46590E-03 0.03470  1.60191E-03 0.05545  5.93843E-04 0.09941 ];
LAMBDA                    (idx, [1:  14]) = [  4.93666E-01 0.03470  1.33611E-02 0.00047  3.25442E-02 0.00081  1.21141E-01 0.00033  3.06818E-01 0.00094  8.67280E-01 0.00150  2.91355E+00 0.00181 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'May 19 2021 21:58:01' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 27])  = 'hexassembly_2D_absorber.inp' ;
WORKING_DIRECTORY         (idx, [1: 26])  = '/nfs/home/jgj/projects/UPV' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode01' ;
CPU_TYPE                  (idx, [1: 46])  = 'AMD Ryzen Threadripper 3970X 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 137367609.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct  4 21:34:37 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct  4 21:35:17 2022' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1664912077112 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 30 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.02283E+00  1.02800E+00  1.02219E+00  1.02498E+00  1.04867E+00  1.05263E+00  1.02384E+00  7.91773E-01  7.91122E-01  1.02798E+00  1.05287E+00  1.02005E+00  7.89558E-01  1.02426E+00  1.05196E+00  1.03010E+00  1.01629E+00  1.04779E+00  1.02655E+00  1.05345E+00  7.90504E-01  1.05409E+00  1.02574E+00  1.01319E+00  1.01636E+00  1.02654E+00  1.02312E+00  1.02609E+00  1.02508E+00  1.05239E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 35])  = '/opt/xsdata/xsdir_files/xsdata_endf' ;
DECAY_DATA_FILE_PATH      (idx, [1: 48])  = '/opt/xsdata/JEFF33/Decay_Data/JEFF33-rdd_all.asc' ;
SFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
NFY_DATA_FILE_PATH        (idx, [1: 48])  = '/opt/xsdata/JEFF33/Fission_Yields/JEFF33-nfy.asc' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 3.7E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.10661E-01 0.00051  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  8.89339E-01 6.4E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.79010E-01 0.00014  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  5.25298E-01 0.00012  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.78928E+00 0.00027  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  1.68602E+01 0.00035  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  1.60739E+01 0.00037  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.45258E+01 0.00051  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.36088E+00 0.00035  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 9997390 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99739E+04 0.00130 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.49683E+01 ;
RUNNING_TIME              (idx, 1)        =  6.65350E-01 ;
INIT_TIME                 (idx, [1:  2])  = [  4.35167E-02  4.35167E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.68333E-03  1.68333E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  6.19433E-01  6.19433E-01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  6.64300E-01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 22.49693 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.41043E+01 0.00057 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  7.00759E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128701.35 ;
ALLOC_MEMSIZE             (idx, 1)        = 1400.32;
MEMSIZE                   (idx, 1)        = 1091.76;
XS_MEMSIZE                (idx, 1)        = 387.35;
MAT_MEMSIZE               (idx, 1)        = 27.72;
RES_MEMSIZE               (idx, 1)        = 8.37;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 668.32;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 308.56;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 23 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 241498 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 2 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 14 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 14 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 309 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.99831E-06 0.00013  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.62158E-01 0.00173 ];
U235_FISS                 (idx, [1:   4]) = [  9.00204E-02 0.00098  9.12290E-01 0.00029 ];
U238_FISS                 (idx, [1:   4]) = [  8.65552E-03 0.00345  8.77097E-02 0.00305 ];
U235_CAPT                 (idx, [1:   4]) = [  2.86726E-02 0.00179  2.49119E-01 0.00155 ];
U238_CAPT                 (idx, [1:   4]) = [  4.29957E-02 0.00150  3.73560E-01 0.00112 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 9997390 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 3.65920E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 1150434 1.15117E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 986249 9.86921E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 7860707 7.86557E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 9997390 1.00037E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -8.17701E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  3.20185E-12 0.00076 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  2.44351E-01 0.00075 ];
TOT_FISSRATE              (idx, [1:   2]) = [  9.86052E-02 0.00076 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  1.14971E-01 0.00061 ];
TOT_ABSRATE               (idx, [1:   2]) = [  2.13577E-01 0.00060 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.99831E-01 0.00013 ];
TOT_FLUX                  (idx, [1:   2]) = [  3.83988E+01 0.00029 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.86423E-01 0.00016 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  1.60777E+01 0.00034 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.86402E+00 0.00073 ];
SIX_FF_F                  (idx, [1:   2]) = [  8.52868E-01 0.00049 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.82670E-01 0.00110 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.55022E+00 0.00120 ];
SIX_FF_LF                 (idx, [1:   2]) = [  2.80565E-01 0.00059 ];
SIX_FF_LT                 (idx, [1:   2]) = [  7.60762E-01 0.00034 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14585E+00 0.00080 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  2.44572E-01 0.00097 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.47807E+00 1.6E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02671E+02 1.7E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  2.44561E-01 0.00098  2.42443E-01 0.00097  2.12865E-03 0.01050 ];
IMP_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
COL_KEFF                  (idx, [1:   2]) = [  2.44393E-01 0.00080 ];
ABS_KEFF                  (idx, [1:   2]) = [  2.44437E-01 0.00075 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14599E+00 0.00034 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.31676E+01 0.00048 ];
IMP_ALF                   (idx, [1:   2]) = [  1.31807E+01 0.00027 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  3.83063E-05 0.00631 ];
IMP_EALF                  (idx, [1:   2]) = [  3.77563E-05 0.00355 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  3.83173E-01 0.00272 ];
IMP_AFGE                  (idx, [1:   2]) = [  3.82479E-01 0.00091 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  3.01927E-02 0.00419  9.24995E-04 0.02176  5.03528E-03 0.00845  4.99499E-03 0.00923  1.16770E-02 0.00547  5.34205E-03 0.00867  2.21837E-03 0.01298 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  5.11675E-01 0.00488  1.33611E-02 0.00014  3.25403E-02 0.00013  1.21199E-01 7.0E-05  3.07058E-01 0.00018  8.66487E-01 0.00026  2.90811E+00 0.00043 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  8.81649E-03 0.00816  2.72268E-04 0.05149  1.50312E-03 0.01995  1.45803E-03 0.02070  3.35929E-03 0.01468  1.58926E-03 0.01893  6.34525E-04 0.03381 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  5.08559E-01 0.01310  1.33544E-02 0.00022  3.25614E-02 0.00032  1.21174E-01 0.00017  3.06933E-01 0.00041  8.66936E-01 0.00062  2.90992E+00 0.00097 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  4.85774E-05 0.00236  4.85134E-05 0.00235  5.61001E-05 0.02312 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.18788E-05 0.00207  1.18632E-05 0.00207  1.37153E-05 0.02293 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  8.68974E-03 0.01043  2.70582E-04 0.06334  1.47407E-03 0.02434  1.47764E-03 0.02685  3.28101E-03 0.01582  1.55596E-03 0.02523  6.30468E-04 0.04396 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  5.07907E-01 0.01635  1.33534E-02 0.00034  3.25553E-02 0.00040  1.21201E-01 0.00021  3.07118E-01 0.00056  8.67255E-01 0.00084  2.90641E+00 0.00138 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.85141E-05 0.00658  4.84099E-05 0.00652  5.75531E-05 0.05433 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.18632E-05 0.00647  1.18378E-05 0.00641  1.40664E-05 0.05415 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  8.66956E-03 0.04306  4.20749E-04 0.22546  1.41031E-03 0.09430  1.20524E-03 0.09711  3.31507E-03 0.07513  1.38061E-03 0.11194  9.37594E-04 0.14050 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  5.79172E-01 0.06119  1.33570E-02 0.00109  3.25956E-02 0.00116  1.21224E-01 0.00076  3.06904E-01 0.00173  8.66676E-01 0.00283  2.90741E+00 0.00357 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  8.72988E-03 0.04178  4.15360E-04 0.22586  1.38787E-03 0.09259  1.24556E-03 0.09633  3.32720E-03 0.07154  1.40820E-03 0.11097  9.45694E-04 0.13655 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  5.86308E-01 0.06010  1.33570E-02 0.00109  3.25985E-02 0.00115  1.21217E-01 0.00074  3.06821E-01 0.00170  8.66841E-01 0.00284  2.90677E+00 0.00356 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.80255E+02 0.04373 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.84821E-05 0.00148 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.18556E-05 0.00104 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  8.55066E-03 0.00690 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.76396E+02 0.00696 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  8.86765E-08 0.00092 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.72310E-06 0.00091  3.72311E-06 0.00093  3.71805E-06 0.00789 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.86564E-05 0.00093  1.86590E-05 0.00094  1.84121E-05 0.00995 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  1.28519E-01 0.00095  1.31202E-01 0.00095  4.24282E-02 0.00932 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.01810E+01 0.00891 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  1.60739E+01 0.00037  2.57880E+01 0.00064 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '7' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  6.37662E+04 0.00831  3.09193E+05 0.00388  7.05072E+05 0.00379  9.65658E+05 0.00413  9.46522E+05 0.00374  9.37846E+05 0.00390  5.59377E+05 0.00316  4.75560E+05 0.00203  5.55946E+05 0.00267  4.42339E+05 0.00265  3.63069E+05 0.00280  2.94869E+05 0.00346  2.60403E+05 0.00426  2.13862E+05 0.00293  1.95961E+05 0.00476  1.35902E+05 0.00281  6.80747E+04 0.00356  2.10870E+05 0.00243  1.72281E+05 0.00248  2.91738E+05 0.00198  2.58449E+05 0.00086  1.68473E+05 0.00305  9.22841E+04 0.00490  9.55905E+04 0.00181  8.21893E+04 0.00590  6.82090E+04 0.00242  1.02943E+05 0.00268  2.22280E+04 0.00309  2.69820E+04 0.00363  2.43394E+04 0.00570  1.35986E+04 0.00725  2.34096E+04 0.01053  1.52012E+04 0.00864  1.19455E+04 0.01367  2.18650E+03 0.01195  2.18895E+03 0.02087  2.22287E+03 0.00915  2.28421E+03 0.00910  2.24995E+03 0.01269  2.20131E+03 0.02352  2.26997E+03 0.02035  2.13906E+03 0.01098  4.05664E+03 0.01459  6.27388E+03 0.01366  7.88496E+03 0.00765  2.04304E+04 0.00963  2.10697E+04 0.00876  2.21933E+04 0.00614  1.41372E+04 0.00857  9.87595E+03 0.01185  7.27919E+03 0.01433  7.98262E+03 0.01152  1.35908E+04 0.00828  1.60191E+04 0.00578  2.60309E+04 0.00412  3.12373E+04 0.00437  3.54736E+04 0.00230  1.85052E+04 0.00485  1.16523E+04 0.00651  7.78121E+03 0.00915  6.50655E+03 0.00622  6.03899E+03 0.00253  4.80222E+03 0.00770  3.10654E+03 0.00743  2.76748E+03 0.01068  2.41390E+03 0.01387  1.96899E+03 0.00643  1.49046E+03 0.01568  9.53848E+02 0.00680  3.39689E+02 0.02635 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.19641E+00 0.00142 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.61449E+00 0.00284  1.36586E-01 0.00396 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.06594E-01 0.00038  7.40894E-01 0.00084 ];
INF_CAPT                  (idx, [1:   4]) = [  2.25383E-03 0.00073  1.68849E-02 0.00081 ];
INF_ABS                   (idx, [1:   4]) = [  3.84314E-03 0.00054  4.97100E-02 0.00097 ];
INF_FISS                  (idx, [1:   4]) = [  1.58930E-03 0.00117  3.28250E-02 0.00117 ];
INF_NSF                   (idx, [1:   4]) = [  3.98252E-03 0.00113  7.99848E-02 0.00117 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.50583E+00 7.7E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02940E+02 6.1E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  2.65491E-08 0.00084  1.96064E-06 0.00148 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  4.02744E-01 0.00039  6.91292E-01 0.00086 ];
INF_SCATT1                (idx, [1:   4]) = [  1.17710E-01 0.00103  2.04077E-01 0.00111 ];
INF_SCATT2                (idx, [1:   4]) = [  4.64378E-02 0.00108  7.48976E-02 0.00987 ];
INF_SCATT3                (idx, [1:   4]) = [  4.17852E-03 0.01354  2.95626E-02 0.01391 ];
INF_SCATT4                (idx, [1:   4]) = [ -3.83634E-03 0.00905  1.35075E-02 0.02690 ];
INF_SCATT5                (idx, [1:   4]) = [  1.48981E-04 0.16759  7.19463E-03 0.02111 ];
INF_SCATT6                (idx, [1:   4]) = [  2.16717E-03 0.01982  4.41137E-03 0.00611 ];
INF_SCATT7                (idx, [1:   4]) = [  3.07813E-04 0.11819  3.07672E-03 0.04179 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  4.02755E-01 0.00039  6.91292E-01 0.00086 ];
INF_SCATTP1               (idx, [1:   4]) = [  1.17711E-01 0.00103  2.04077E-01 0.00111 ];
INF_SCATTP2               (idx, [1:   4]) = [  4.64375E-02 0.00108  7.48976E-02 0.00987 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.17864E-03 0.01356  2.95626E-02 0.01391 ];
INF_SCATTP4               (idx, [1:   4]) = [ -3.83594E-03 0.00901  1.35075E-02 0.02690 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.48890E-04 0.16836  7.19463E-03 0.02111 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.16697E-03 0.01981  4.41137E-03 0.00611 ];
INF_SCATTP7               (idx, [1:   4]) = [  3.08008E-04 0.11807  3.07672E-03 0.04179 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.32812E-01 0.00043  5.00451E-01 0.00110 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.43177E+00 0.00043  6.66069E-01 0.00110 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  3.83284E-03 0.00052  4.97100E-02 0.00097 ];
INF_REMXS                 (idx, [1:   4]) = [  7.23242E-03 0.00035  5.43079E-02 0.00170 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37701E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49395E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.94297E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.92615E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55062E-04 ;
I135_BR                   (idx, 1)        =  8.34918E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  3.99362E-01 0.00039  3.38272E-03 0.00250  4.70591E-03 0.01634  6.86586E-01 0.00098 ];
INF_S1                    (idx, [1:   8]) = [  1.16699E-01 0.00106  1.01111E-03 0.00390  1.59808E-03 0.02544  2.02479E-01 0.00103 ];
INF_S2                    (idx, [1:   8]) = [  4.66649E-02 0.00111 -2.27061E-04 0.01590  7.78419E-04 0.02291  7.41192E-02 0.00995 ];
INF_S3                    (idx, [1:   8]) = [  4.54074E-03 0.01215 -3.62223E-04 0.00826  2.06149E-04 0.07026  2.93565E-02 0.01396 ];
INF_S4                    (idx, [1:   8]) = [ -3.69412E-03 0.00956 -1.42221E-04 0.02105 -1.49076E-05 1.00000  1.35224E-02 0.02754 ];
INF_S5                    (idx, [1:   8]) = [  1.71430E-04 0.14451 -2.24496E-05 0.09939 -1.21022E-04 0.13047  7.31565E-03 0.01983 ];
INF_S6                    (idx, [1:   8]) = [  2.17295E-03 0.02000 -5.78194E-06 0.56863 -6.92908E-05 0.12784  4.48066E-03 0.00627 ];
INF_S7                    (idx, [1:   8]) = [  3.12820E-04 0.11521 -5.00678E-06 0.54899 -4.34181E-05 0.34890  3.12014E-03 0.03965 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.99372E-01 0.00039  3.38272E-03 0.00250  4.70591E-03 0.01634  6.86586E-01 0.00098 ];
INF_SP1                   (idx, [1:   8]) = [  1.16699E-01 0.00106  1.01111E-03 0.00390  1.59808E-03 0.02544  2.02479E-01 0.00103 ];
INF_SP2                   (idx, [1:   8]) = [  4.66646E-02 0.00111 -2.27061E-04 0.01590  7.78419E-04 0.02291  7.41192E-02 0.00995 ];
INF_SP3                   (idx, [1:   8]) = [  4.54087E-03 0.01217 -3.62223E-04 0.00826  2.06149E-04 0.07026  2.93565E-02 0.01396 ];
INF_SP4                   (idx, [1:   8]) = [ -3.69372E-03 0.00952 -1.42221E-04 0.02105 -1.49076E-05 1.00000  1.35224E-02 0.02754 ];
INF_SP5                   (idx, [1:   8]) = [  1.71340E-04 0.14514 -2.24496E-05 0.09939 -1.21022E-04 0.13047  7.31565E-03 0.01983 ];
INF_SP6                   (idx, [1:   8]) = [  2.17275E-03 0.01999 -5.78194E-06 0.56863 -6.92908E-05 0.12784  4.48066E-03 0.00627 ];
INF_SP7                   (idx, [1:   8]) = [  3.13015E-04 0.11511 -5.00678E-06 0.54899 -4.34181E-05 0.34890  3.12014E-03 0.03965 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  4.43343E-02 0.00255  8.01635E-02 0.00369 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  5.08484E-02 0.00241  1.50824E-01 0.00924 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  5.02170E-02 0.00288  7.37278E-02 0.00574 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  3.56028E-02 0.00255  5.80511E-02 0.00496 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  7.51882E+00 0.00255  4.15839E+00 0.00368 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  6.55559E+00 0.00241  2.21083E+00 0.00911 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  6.63808E+00 0.00287  4.52173E+00 0.00572 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  9.36279E+00 0.00255  5.74263E+00 0.00496 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  8.92541E-03 0.02585  2.96634E-04 0.13524  1.45232E-03 0.06739  1.49859E-03 0.05748  3.43221E-03 0.04412  1.58379E-03 0.06359  6.61871E-04 0.07952 ];
LAMBDA                    (idx, [1:  14]) = [  5.12090E-01 0.03149  1.33593E-02 0.00050  3.25815E-02 0.00061  1.21090E-01 0.00033  3.06532E-01 0.00088  8.70871E-01 0.00148  2.90332E+00 0.00195 ];
