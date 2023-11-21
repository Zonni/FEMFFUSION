
% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50812E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19149 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95769E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '1' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.21617E+06 0.00175  6.16802E+06 0.00088  1.44108E+07 0.00093  2.18694E+07 0.00058  2.51174E+07 0.00091  2.72449E+07 0.00072  1.55633E+07 0.00077  1.35497E+07 0.00074  2.35691E+07 0.00068  2.19213E+07 0.00058  2.43825E+07 0.00076  2.33413E+07 0.00077  2.31357E+07 0.00064  1.68338E+07 0.00079  1.26557E+07 0.00105  6.40663E+06 0.00106  1.62785E+06 0.00098  6.93502E+06 0.00085  6.55130E+06 0.00086  7.57657E+06 0.00093  3.69417E+06 0.00074  1.69295E+06 0.00117  8.43653E+05 0.00184  7.52809E+05 0.00212  6.50681E+05 0.00233  6.21233E+05 0.00153  9.82954E+05 0.00196  3.01948E+05 0.00298  3.94246E+05 0.00255  4.17785E+05 0.00172  2.26068E+05 0.00325  4.31694E+05 0.00214  2.86083E+05 0.00262  1.90269E+05 0.00154  2.82828E+04 0.00191  2.72450E+04 0.00414  2.83954E+04 0.00316  2.92131E+04 0.00457  2.93033E+04 0.00570  2.94626E+04 0.00373  3.02657E+04 0.00245  2.88891E+04 0.00410  5.42539E+04 0.00267  8.66204E+04 0.00285  1.10667E+05 0.00288  2.98456E+05 0.00160  3.46892E+05 0.00101  4.18134E+05 0.00083  2.78705E+05 0.00069  1.93221E+05 0.00153  1.43178E+05 0.00241  1.59244E+05 0.00187  2.84961E+05 0.00033  3.59363E+05 0.00105  6.15685E+05 0.00158  7.43577E+05 0.00062  7.82027E+05 0.00085  3.58366E+05 0.00094  2.05959E+05 0.00068  1.25116E+05 0.00218  9.72204E+04 0.00307  8.23936E+04 0.00119  5.90653E+04 0.00265  3.44424E+04 0.00156  2.80241E+04 0.00223  2.17376E+04 0.00222  1.53736E+04 0.00447  9.57718E+03 0.00214  4.60778E+03 0.00197  1.02502E+03 0.00900 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.54256E+00 0.00032 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  1.17459E+18 0.00077  2.01862E+16 0.00051 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  3.50422E-01 2.1E-05  5.74658E-01 8.4E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  3.64125E-03 0.00029  5.71914E-02 0.00013 ];
INF_ABS                   (idx, [1:   4]) = [  7.44673E-03 0.00029  2.90015E-01 0.00016 ];
INF_FISS                  (idx, [1:   4]) = [  3.80547E-03 0.00030  2.32824E-01 0.00017 ];
INF_NSF                   (idx, [1:   4]) = [  9.41276E-03 0.00029  5.67321E-01 0.00017 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.47348E+00 1.1E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02531E+02 7.4E-07  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  1.28636E-08 0.00060  1.82870E-06 0.00013 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.42975E-01 2.6E-05  2.84652E-01 0.00029 ];
INF_SCATT1                (idx, [1:   4]) = [  3.09440E-02 9.5E-05  7.77698E-03 0.00643 ];
INF_SCATT2                (idx, [1:   4]) = [  1.13849E-02 0.00033  2.55328E-04 0.26745 ];
INF_SCATT3                (idx, [1:   4]) = [  3.18270E-03 0.00077  4.25132E-05 0.98355 ];
INF_SCATT4                (idx, [1:   4]) = [  1.58550E-03 0.00254 -4.96218E-05 0.97832 ];
INF_SCATT5                (idx, [1:   4]) = [  5.31825E-04 0.00705 -5.83736E-05 0.25265 ];
INF_SCATT6                (idx, [1:   4]) = [  2.35883E-04 0.01633  2.60417E-06 1.00000 ];
INF_SCATT7                (idx, [1:   4]) = [  8.93911E-05 0.08355  2.52555E-05 0.49105 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.42984E-01 2.6E-05  2.84652E-01 0.00029 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.09441E-02 9.4E-05  7.77698E-03 0.00643 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.13849E-02 0.00033  2.55328E-04 0.26745 ];
INF_SCATTP3               (idx, [1:   4]) = [  3.18269E-03 0.00077  4.25132E-05 0.98355 ];
INF_SCATTP4               (idx, [1:   4]) = [  1.58546E-03 0.00255 -4.96218E-05 0.97832 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.31816E-04 0.00706 -5.83736E-05 0.25265 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.35879E-04 0.01638  2.60417E-06 1.00000 ];
INF_SCATTP7               (idx, [1:   4]) = [  8.94032E-05 0.08360  2.52555E-05 0.49105 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.82711E-01 3.2E-05  5.40531E-01 0.00013 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.17906E+00 3.2E-05  6.16678E-01 0.00013 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  7.43753E-03 0.00030  2.90015E-01 0.00016 ];
INF_REMXS                 (idx, [1:   4]) = [  7.55589E-03 0.00023  2.92173E-01 0.00033 ];

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

INF_S0                    (idx, [1:   8]) = [  3.42866E-01 2.6E-05  1.08087E-04 0.00092  2.16641E-03 0.00501  2.82485E-01 0.00028 ];
INF_S1                    (idx, [1:   8]) = [  3.09693E-02 9.4E-05 -2.52865E-05 0.00520 -1.87516E-04 0.01844  7.96450E-03 0.00618 ];
INF_S2                    (idx, [1:   8]) = [  1.13875E-02 0.00032 -2.63191E-06 0.04348 -1.02925E-04 0.02414  3.58252E-04 0.18634 ];
INF_S3                    (idx, [1:   8]) = [  3.18350E-03 0.00076 -8.01423E-07 0.07570 -4.05382E-05 0.03976  8.30514E-05 0.51197 ];
INF_S4                    (idx, [1:   8]) = [  1.58585E-03 0.00257 -3.52464E-07 0.14679 -1.52107E-05 0.17231 -3.44111E-05 1.00000 ];
INF_S5                    (idx, [1:   8]) = [  5.31956E-04 0.00711 -1.30554E-07 0.83030 -8.63408E-06 0.20483 -4.97395E-05 0.32329 ];
INF_S6                    (idx, [1:   8]) = [  2.35967E-04 0.01668 -8.36449E-08 1.00000 -7.75402E-06 0.36323  1.03582E-05 1.00000 ];
INF_S7                    (idx, [1:   8]) = [  8.94218E-05 0.08328 -3.06609E-08 1.00000 -6.51059E-06 0.39022  3.17661E-05 0.41252 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.42876E-01 2.6E-05  1.08087E-04 0.00092  2.16641E-03 0.00501  2.82485E-01 0.00028 ];
INF_SP1                   (idx, [1:   8]) = [  3.09694E-02 9.3E-05 -2.52865E-05 0.00520 -1.87516E-04 0.01844  7.96450E-03 0.00618 ];
INF_SP2                   (idx, [1:   8]) = [  1.13875E-02 0.00032 -2.63191E-06 0.04348 -1.02925E-04 0.02414  3.58252E-04 0.18634 ];
INF_SP3                   (idx, [1:   8]) = [  3.18349E-03 0.00076 -8.01423E-07 0.07570 -4.05382E-05 0.03976  8.30514E-05 0.51197 ];
INF_SP4                   (idx, [1:   8]) = [  1.58582E-03 0.00257 -3.52464E-07 0.14679 -1.52107E-05 0.17231 -3.44111E-05 1.00000 ];
INF_SP5                   (idx, [1:   8]) = [  5.31946E-04 0.00712 -1.30554E-07 0.83030 -8.63408E-06 0.20483 -4.97395E-05 0.32329 ];
INF_SP6                   (idx, [1:   8]) = [  2.35963E-04 0.01672 -8.36449E-08 1.00000 -7.75402E-06 0.36323  1.03582E-05 1.00000 ];
INF_SP7                   (idx, [1:   8]) = [  8.94339E-05 0.08333 -3.06609E-08 1.00000 -6.51059E-06 0.39022  3.17661E-05 0.41252 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  3.95138E-01 0.00048  6.33600E-01 0.00523 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  4.22326E-01 0.00054  1.02169E+00 0.01031 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  4.22479E-01 0.00059  1.03044E+00 0.01265 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  3.49961E-01 0.00068  3.59070E-01 0.00335 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  8.43588E-01 0.00048  5.26152E-01 0.00526 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  7.89280E-01 0.00054  3.26396E-01 0.01034 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  7.88995E-01 0.00059  3.23695E-01 0.01278 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  9.52489E-01 0.00068  9.28365E-01 0.00335 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  7.10181E-03 0.00332  2.40430E-04 0.01834  1.23765E-03 0.00891  1.21059E-03 0.00743  2.76071E-03 0.00533  1.16484E-03 0.00838  4.87594E-04 0.01323 ];
LAMBDA                    (idx, [1:  14]) = [  4.82586E-01 0.00488  1.33444E-02 6.4E-05  3.26760E-02 6.8E-05  1.20916E-01 3.8E-05  3.04259E-01 9.3E-05  8.55691E-01 0.00015  2.87322E+00 0.00027 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50812E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19148 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95769E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '2' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  4.89844E+05 0.00190  2.48794E+06 0.00259  5.79933E+06 0.00151  8.81537E+06 0.00182  1.01237E+07 0.00166  1.09972E+07 0.00146  6.28746E+06 0.00137  5.47963E+06 0.00098  9.58043E+06 0.00164  8.94285E+06 0.00160  1.00113E+07 0.00176  9.66232E+06 0.00184  9.68006E+06 0.00185  7.13224E+06 0.00139  5.42869E+06 0.00181  2.76821E+06 0.00164  7.02932E+05 0.00189  3.00299E+06 0.00169  2.86390E+06 0.00174  3.31937E+06 0.00182  1.57330E+06 0.00274  6.96068E+05 0.00150  3.42529E+05 0.00320  3.02667E+05 0.00320  2.59866E+05 0.00256  2.48065E+05 0.00570  3.93219E+05 0.00398  1.19946E+05 0.00487  1.57160E+05 0.00244  1.65251E+05 0.00221  8.99395E+04 0.00428  1.71281E+05 0.00428  1.13348E+05 0.00327  7.60807E+04 0.00179  1.12032E+04 0.00919  1.08895E+04 0.00993  1.11900E+04 0.00254  1.16418E+04 0.00737  1.16793E+04 0.00920  1.15328E+04 0.00301  1.20428E+04 0.00928  1.13102E+04 0.00810  2.13583E+04 0.00372  3.43041E+04 0.00367  4.36542E+04 0.00470  1.17998E+05 0.00331  1.33783E+05 0.00378  1.59154E+05 0.00181  1.05495E+05 0.00094  7.32958E+04 0.00168  5.37474E+04 0.00103  5.97289E+04 0.00441  1.06598E+05 0.00065  1.33950E+05 0.00232  2.28753E+05 0.00185  2.76872E+05 0.00217  2.91340E+05 0.00201  1.33077E+05 0.00260  7.62898E+04 0.00124  4.64512E+04 0.00306  3.60696E+04 0.00316  3.05691E+04 0.00190  2.18783E+04 0.00359  1.29005E+04 0.00328  1.04941E+04 0.00405  8.05161E+03 0.00616  5.73309E+03 0.00495  3.58646E+03 0.00655  1.69445E+03 0.00590  3.82822E+02 0.00535 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.56275E+00 0.00033 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.83581E+17 0.00145  7.55831E+15 0.00137 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  3.50802E-01 8.7E-05  5.73699E-01 7.7E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  3.38042E-03 0.00022  5.70055E-02 0.00010 ];
INF_ABS                   (idx, [1:   4]) = [  7.21325E-03 0.00013  2.89081E-01 0.00015 ];
INF_FISS                  (idx, [1:   4]) = [  3.83283E-03 0.00022  2.32075E-01 0.00016 ];
INF_NSF                   (idx, [1:   4]) = [  9.47614E-03 0.00022  5.65497E-01 0.00016 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.47236E+00 1.4E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02524E+02 9.4E-07  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  1.27458E-08 0.00064  1.82365E-06 0.00016 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.43593E-01 9.1E-05  2.84548E-01 0.00041 ];
INF_SCATT1                (idx, [1:   4]) = [  3.05781E-02 0.00066  7.95270E-03 0.01079 ];
INF_SCATT2                (idx, [1:   4]) = [  1.11770E-02 0.00119  4.21591E-04 0.10389 ];
INF_SCATT3                (idx, [1:   4]) = [  3.11904E-03 0.00192 -4.83243E-06 1.00000 ];
INF_SCATT4                (idx, [1:   4]) = [  1.54346E-03 0.00158 -5.80785E-05 1.00000 ];
INF_SCATT5                (idx, [1:   4]) = [  5.32420E-04 0.02075 -1.06085E-04 0.39023 ];
INF_SCATT6                (idx, [1:   4]) = [  2.35965E-04 0.02262 -3.28752E-05 1.00000 ];
INF_SCATT7                (idx, [1:   4]) = [  8.51368E-05 0.08220 -8.43531E-05 0.34339 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.43602E-01 9.1E-05  2.84548E-01 0.00041 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.05782E-02 0.00066  7.95270E-03 0.01079 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.11769E-02 0.00119  4.21591E-04 0.10389 ];
INF_SCATTP3               (idx, [1:   4]) = [  3.11903E-03 0.00190 -4.83243E-06 1.00000 ];
INF_SCATTP4               (idx, [1:   4]) = [  1.54342E-03 0.00159 -5.80785E-05 1.00000 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.32438E-04 0.02075 -1.06085E-04 0.39023 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.36022E-04 0.02268 -3.28752E-05 1.00000 ];
INF_SCATTP7               (idx, [1:   4]) = [  8.51383E-05 0.08244 -8.43531E-05 0.34339 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.83496E-01 0.00015  5.39328E-01 0.00020 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.17580E+00 0.00015  6.18053E-01 0.00020 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  7.20429E-03 0.00013  2.89081E-01 0.00015 ];
INF_REMXS                 (idx, [1:   4]) = [  7.31291E-03 0.00016  2.91410E-01 0.00040 ];

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

INF_S0                    (idx, [1:   8]) = [  3.43490E-01 9.1E-05  1.03786E-04 0.00286  2.25881E-03 0.00980  2.82289E-01 0.00035 ];
INF_S1                    (idx, [1:   8]) = [  3.06023E-02 0.00066 -2.41825E-05 0.00645 -1.95116E-04 0.06341  8.14781E-03 0.01071 ];
INF_S2                    (idx, [1:   8]) = [  1.11798E-02 0.00120 -2.82254E-06 0.04438 -1.07775E-04 0.04491  5.29365E-04 0.08895 ];
INF_S3                    (idx, [1:   8]) = [  3.11981E-03 0.00191 -7.74988E-07 0.20125 -4.43094E-05 0.19654  3.94770E-05 1.00000 ];
INF_S4                    (idx, [1:   8]) = [  1.54335E-03 0.00161  1.09902E-07 1.00000 -1.78138E-05 0.13121 -4.02647E-05 1.00000 ];
INF_S5                    (idx, [1:   8]) = [  5.32756E-04 0.02062 -3.36793E-07 0.58066 -1.06008E-05 0.23968 -9.54847E-05 0.43410 ];
INF_S6                    (idx, [1:   8]) = [  2.36070E-04 0.02283 -1.04527E-07 0.53565 -5.68501E-07 1.00000 -3.23067E-05 1.00000 ];
INF_S7                    (idx, [1:   8]) = [  8.51426E-05 0.08143 -5.79998E-09 1.00000 -3.06598E-06 1.00000 -8.12871E-05 0.31937 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.43499E-01 9.1E-05  1.03786E-04 0.00286  2.25881E-03 0.00980  2.82289E-01 0.00035 ];
INF_SP1                   (idx, [1:   8]) = [  3.06024E-02 0.00066 -2.41825E-05 0.00645 -1.95116E-04 0.06341  8.14781E-03 0.01071 ];
INF_SP2                   (idx, [1:   8]) = [  1.11797E-02 0.00120 -2.82254E-06 0.04438 -1.07775E-04 0.04491  5.29365E-04 0.08895 ];
INF_SP3                   (idx, [1:   8]) = [  3.11981E-03 0.00189 -7.74988E-07 0.20125 -4.43094E-05 0.19654  3.94770E-05 1.00000 ];
INF_SP4                   (idx, [1:   8]) = [  1.54331E-03 0.00162  1.09902E-07 1.00000 -1.78138E-05 0.13121 -4.02647E-05 1.00000 ];
INF_SP5                   (idx, [1:   8]) = [  5.32774E-04 0.02062 -3.36793E-07 0.58066 -1.06008E-05 0.23968 -9.54847E-05 0.43410 ];
INF_SP6                   (idx, [1:   8]) = [  2.36126E-04 0.02288 -1.04527E-07 0.53565 -5.68501E-07 1.00000 -3.23067E-05 1.00000 ];
INF_SP7                   (idx, [1:   8]) = [  8.51441E-05 0.08167 -5.79998E-09 1.00000 -3.06598E-06 1.00000 -8.12871E-05 0.31937 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  1.13620E-01 0.00122  9.16665E-02 0.00199 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  1.17843E-01 0.00130  6.34249E-02 0.00228 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  1.23367E-01 0.00128  1.48366E-01 0.00441 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  1.01915E-01 0.00114  9.78466E-02 0.00191 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  2.93378E+00 0.00122  3.63643E+00 0.00199 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  2.82864E+00 0.00130  5.25567E+00 0.00228 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  2.70199E+00 0.00128  2.24687E+00 0.00441 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  3.27071E+00 0.00114  3.40674E+00 0.00192 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  7.04632E-03 0.00508  2.32010E-04 0.02915  1.25154E-03 0.01348  1.19527E-03 0.01266  2.70783E-03 0.00821  1.17994E-03 0.01255  4.79730E-04 0.02156 ];
LAMBDA                    (idx, [1:  14]) = [  4.82785E-01 0.00838  1.33460E-02 0.00014  3.26781E-02 0.00010  1.20925E-01 6.4E-05  3.04255E-01 0.00013  8.56037E-01 0.00030  2.87355E+00 0.00042 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50812E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19147 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95768E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '3' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  4.89999E+05 0.00176  2.48113E+06 0.00295  5.77981E+06 0.00175  8.78681E+06 0.00207  1.00956E+07 0.00185  1.09757E+07 0.00179  6.27518E+06 0.00169  5.46927E+06 0.00171  9.56165E+06 0.00143  8.92094E+06 0.00170  9.99176E+06 0.00157  9.63920E+06 0.00158  9.66531E+06 0.00203  7.12532E+06 0.00194  5.42407E+06 0.00239  2.76290E+06 0.00233  7.01470E+05 0.00239  2.99628E+06 0.00207  2.85542E+06 0.00188  3.30866E+06 0.00195  1.57282E+06 0.00263  6.95161E+05 0.00223  3.39609E+05 0.00307  3.00081E+05 0.00326  2.60136E+05 0.00429  2.49571E+05 0.00452  3.89273E+05 0.00213  1.19236E+05 0.00275  1.55993E+05 0.00274  1.65975E+05 0.00168  8.98996E+04 0.00272  1.71657E+05 0.00314  1.13176E+05 0.00266  7.59808E+04 0.00361  1.12711E+04 0.00500  1.09214E+04 0.00720  1.11758E+04 0.00339  1.16704E+04 0.00848  1.16041E+04 0.00546  1.14638E+04 0.01040  1.20061E+04 0.00849  1.14464E+04 0.00867  2.14989E+04 0.00336  3.44454E+04 0.00574  4.41617E+04 0.00604  1.17988E+05 0.00561  1.34094E+05 0.00139  1.58846E+05 0.00149  1.05486E+05 0.00301  7.28397E+04 0.00524  5.37170E+04 0.00160  5.96401E+04 0.00270  1.06525E+05 0.00158  1.34154E+05 0.00319  2.29008E+05 0.00249  2.76592E+05 0.00269  2.90695E+05 0.00200  1.32690E+05 0.00225  7.63452E+04 0.00290  4.62247E+04 0.00303  3.58728E+04 0.00319  3.05905E+04 0.00294  2.19052E+04 0.00553  1.28718E+04 0.00640  1.05053E+04 0.00583  7.97289E+03 0.00313  5.73103E+03 0.00561  3.60247E+03 0.00446  1.69551E+03 0.00474  3.83405E+02 0.00638 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.56247E+00 0.00019 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.82534E+17 0.00171  7.55114E+15 0.00192 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  3.50820E-01 5.6E-05  5.73630E-01 9.2E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  3.38004E-03 0.00042  5.69954E-02 0.00018 ];
INF_ABS                   (idx, [1:   4]) = [  7.21354E-03 0.00042  2.89013E-01 0.00018 ];
INF_FISS                  (idx, [1:   4]) = [  3.83350E-03 0.00044  2.32018E-01 0.00018 ];
INF_NSF                   (idx, [1:   4]) = [  9.47774E-03 0.00043  5.65357E-01 0.00018 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.47235E+00 1.1E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02524E+02 9.4E-07  2.02270E+02 5.9E-09 ];
INF_INVV                  (idx, [1:   4]) = [  1.27539E-08 0.00029  1.82329E-06 0.00014 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.43608E-01 5.6E-05  2.84488E-01 0.00028 ];
INF_SCATT1                (idx, [1:   4]) = [  3.05719E-02 0.00017  7.69150E-03 0.02722 ];
INF_SCATT2                (idx, [1:   4]) = [  1.11701E-02 0.00095  2.56215E-04 0.14842 ];
INF_SCATT3                (idx, [1:   4]) = [  3.12272E-03 0.00295  8.89952E-05 0.85257 ];
INF_SCATT4                (idx, [1:   4]) = [  1.55230E-03 0.00629 -2.72200E-05 1.00000 ];
INF_SCATT5                (idx, [1:   4]) = [  5.17925E-04 0.00592 -2.95364E-05 1.00000 ];
INF_SCATT6                (idx, [1:   4]) = [  2.34638E-04 0.02558 -9.11618E-06 1.00000 ];
INF_SCATT7                (idx, [1:   4]) = [  8.05960E-05 0.04575 -2.12374E-05 1.00000 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.43617E-01 5.6E-05  2.84488E-01 0.00028 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.05721E-02 0.00017  7.69150E-03 0.02722 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.11701E-02 0.00095  2.56215E-04 0.14842 ];
INF_SCATTP3               (idx, [1:   4]) = [  3.12269E-03 0.00295  8.89952E-05 0.85257 ];
INF_SCATTP4               (idx, [1:   4]) = [  1.55236E-03 0.00631 -2.72200E-05 1.00000 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.17940E-04 0.00593 -2.95364E-05 1.00000 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.34633E-04 0.02554 -9.11618E-06 1.00000 ];
INF_SCATTP7               (idx, [1:   4]) = [  8.06170E-05 0.04576 -2.12374E-05 1.00000 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.83550E-01 4.8E-05  5.39510E-01 0.00048 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.17557E+00 4.8E-05  6.17845E-01 0.00048 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  7.20451E-03 0.00043  2.89013E-01 0.00018 ];
INF_REMXS                 (idx, [1:   4]) = [  7.31588E-03 0.00023  2.91390E-01 0.00027 ];

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

INF_S0                    (idx, [1:   8]) = [  3.43504E-01 5.7E-05  1.03413E-04 0.00612  2.24883E-03 0.00607  2.82240E-01 0.00031 ];
INF_S1                    (idx, [1:   8]) = [  3.05963E-02 0.00017 -2.44315E-05 0.00816 -1.94645E-04 0.02812  7.88614E-03 0.02590 ];
INF_S2                    (idx, [1:   8]) = [  1.11726E-02 0.00095 -2.58149E-06 0.06530 -1.09805E-04 0.06246  3.66021E-04 0.09217 ];
INF_S3                    (idx, [1:   8]) = [  3.12328E-03 0.00296 -5.57235E-07 0.18175 -4.35847E-05 0.12474  1.32580E-04 0.58869 ];
INF_S4                    (idx, [1:   8]) = [  1.55255E-03 0.00624 -2.50680E-07 0.37747 -1.21430E-05 0.27174 -1.50770E-05 1.00000 ];
INF_S5                    (idx, [1:   8]) = [  5.18280E-04 0.00589 -3.54716E-07 0.28161 -8.80763E-06 0.33650 -2.07287E-05 1.00000 ];
INF_S6                    (idx, [1:   8]) = [  2.34858E-04 0.02537 -2.20641E-07 0.36296 -2.47806E-06 1.00000 -6.63813E-06 1.00000 ];
INF_S7                    (idx, [1:   8]) = [  8.05123E-05 0.04603  8.36942E-08 1.00000  8.45918E-07 1.00000 -2.20833E-05 1.00000 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.43514E-01 5.7E-05  1.03413E-04 0.00612  2.24883E-03 0.00607  2.82240E-01 0.00031 ];
INF_SP1                   (idx, [1:   8]) = [  3.05965E-02 0.00018 -2.44315E-05 0.00816 -1.94645E-04 0.02812  7.88614E-03 0.02590 ];
INF_SP2                   (idx, [1:   8]) = [  1.11727E-02 0.00095 -2.58149E-06 0.06530 -1.09805E-04 0.06246  3.66021E-04 0.09217 ];
INF_SP3                   (idx, [1:   8]) = [  3.12325E-03 0.00296 -5.57235E-07 0.18175 -4.35847E-05 0.12474  1.32580E-04 0.58869 ];
INF_SP4                   (idx, [1:   8]) = [  1.55261E-03 0.00627 -2.50680E-07 0.37747 -1.21430E-05 0.27174 -1.50770E-05 1.00000 ];
INF_SP5                   (idx, [1:   8]) = [  5.18294E-04 0.00591 -3.54716E-07 0.28161 -8.80763E-06 0.33650 -2.07287E-05 1.00000 ];
INF_SP6                   (idx, [1:   8]) = [  2.34854E-04 0.02534 -2.20641E-07 0.36296 -2.47806E-06 1.00000 -6.63813E-06 1.00000 ];
INF_SP7                   (idx, [1:   8]) = [  8.05333E-05 0.04605  8.36942E-08 1.00000  8.45918E-07 1.00000 -2.20833E-05 1.00000 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  8.71406E-02 0.00156  5.67017E-02 0.00271 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  9.36954E-02 0.00177  6.46898E-02 0.00473 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  9.05190E-02 0.00139  4.09721E-02 0.00230 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  7.86980E-02 0.00156  7.66714E-02 0.00206 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  3.82527E+00 0.00156  5.87889E+00 0.00271 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  3.55767E+00 0.00176  5.15326E+00 0.00472 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  3.68250E+00 0.00139  8.13579E+00 0.00229 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  4.23564E+00 0.00156  4.34763E+00 0.00206 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  7.00699E-03 0.00573  2.29536E-04 0.02995  1.22581E-03 0.01391  1.16166E-03 0.01148  2.73191E-03 0.00985  1.18390E-03 0.01367  4.74173E-04 0.01757 ];
LAMBDA                    (idx, [1:  14]) = [  4.84304E-01 0.00744  1.33433E-02 9.2E-05  3.26674E-02 0.00012  1.20920E-01 6.3E-05  3.04275E-01 0.00015  8.55767E-01 0.00024  2.87330E+00 0.00044 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50812E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19147 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95768E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '4' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  4.93091E+05 0.00290  2.48417E+06 0.00207  5.79585E+06 0.00251  8.81819E+06 0.00224  1.01215E+07 0.00244  1.09906E+07 0.00249  6.27984E+06 0.00249  5.47897E+06 0.00250  9.57494E+06 0.00244  8.93586E+06 0.00258  1.00038E+07 0.00279  9.65188E+06 0.00268  9.67274E+06 0.00264  7.12873E+06 0.00253  5.41716E+06 0.00317  2.76452E+06 0.00311  7.01956E+05 0.00332  3.00087E+06 0.00334  2.86293E+06 0.00287  3.31664E+06 0.00240  1.57362E+06 0.00215  6.97745E+05 0.00354  3.40989E+05 0.00305  3.02631E+05 0.00357  2.60791E+05 0.00183  2.47422E+05 0.00615  3.92787E+05 0.00313  1.21202E+05 0.00476  1.57625E+05 0.00393  1.67721E+05 0.00519  9.01436E+04 0.00553  1.71863E+05 0.00267  1.13095E+05 0.00384  7.61094E+04 0.00370  1.12152E+04 0.00386  1.08913E+04 0.00535  1.12098E+04 0.00826  1.16020E+04 0.01391  1.16471E+04 0.00700  1.15471E+04 0.00394  1.21868E+04 0.01184  1.13621E+04 0.00883  2.14819E+04 0.00586  3.42628E+04 0.00622  4.39963E+04 0.00190  1.18638E+05 0.00473  1.34943E+05 0.00435  1.59776E+05 0.00420  1.05585E+05 0.00370  7.24606E+04 0.00233  5.38193E+04 0.00535  5.99026E+04 0.00206  1.06561E+05 0.00274  1.34759E+05 0.00268  2.29447E+05 0.00220  2.76211E+05 0.00286  2.91191E+05 0.00187  1.33175E+05 0.00219  7.66198E+04 0.00268  4.65071E+04 0.00357  3.62893E+04 0.00339  3.07762E+04 0.00591  2.20091E+04 0.00457  1.28352E+04 0.00540  1.04140E+04 0.00301  8.05850E+03 0.00805  5.71487E+03 0.00090  3.58543E+03 0.00745  1.71861E+03 0.00672  3.85816E+02 0.01069 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.56313E+00 0.00016 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.83298E+17 0.00242  7.56902E+15 0.00230 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  3.50760E-01 8.9E-05  5.73628E-01 0.00012 ];
INF_CAPT                  (idx, [1:   4]) = [  3.38098E-03 0.00046  5.69902E-02 0.00026 ];
INF_ABS                   (idx, [1:   4]) = [  7.21594E-03 0.00036  2.89010E-01 0.00024 ];
INF_FISS                  (idx, [1:   4]) = [  3.83496E-03 0.00028  2.32020E-01 0.00024 ];
INF_NSF                   (idx, [1:   4]) = [  9.48149E-03 0.00028  5.65363E-01 0.00024 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.47238E+00 1.0E-05  2.43670E+00 8.3E-09 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02524E+02 8.0E-07  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  1.27737E-08 0.00053  1.82337E-06 0.00019 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.43545E-01 9.5E-05  2.84760E-01 0.00015 ];
INF_SCATT1                (idx, [1:   4]) = [  3.05791E-02 0.00045  7.94785E-03 0.00939 ];
INF_SCATT2                (idx, [1:   4]) = [  1.11803E-02 0.00100  3.03612E-04 0.22021 ];
INF_SCATT3                (idx, [1:   4]) = [  3.11301E-03 0.00467 -3.77080E-05 1.00000 ];
INF_SCATT4                (idx, [1:   4]) = [  1.54777E-03 0.00380 -1.64383E-05 1.00000 ];
INF_SCATT5                (idx, [1:   4]) = [  5.21644E-04 0.01270  4.43723E-06 1.00000 ];
INF_SCATT6                (idx, [1:   4]) = [  2.29616E-04 0.04344 -6.86446E-06 1.00000 ];
INF_SCATT7                (idx, [1:   4]) = [  7.87141E-05 0.07759  6.09884E-05 0.39785 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.43554E-01 9.5E-05  2.84760E-01 0.00015 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.05793E-02 0.00045  7.94785E-03 0.00939 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.11803E-02 0.00100  3.03612E-04 0.22021 ];
INF_SCATTP3               (idx, [1:   4]) = [  3.11301E-03 0.00467 -3.77080E-05 1.00000 ];
INF_SCATTP4               (idx, [1:   4]) = [  1.54773E-03 0.00380 -1.64383E-05 1.00000 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.21683E-04 0.01269  4.43723E-06 1.00000 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.29651E-04 0.04343 -6.86446E-06 1.00000 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.87198E-05 0.07782  6.09884E-05 0.39785 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.83443E-01 0.00014  5.39154E-01 0.00017 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.17602E+00 0.00014  6.18252E-01 0.00017 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  7.20694E-03 0.00035  2.89010E-01 0.00024 ];
INF_REMXS                 (idx, [1:   4]) = [  7.31882E-03 0.00039  2.91122E-01 0.00030 ];

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

INF_S0                    (idx, [1:   8]) = [  3.43441E-01 9.5E-05  1.03696E-04 0.00344  2.25379E-03 0.00735  2.82506E-01 0.00016 ];
INF_S1                    (idx, [1:   8]) = [  3.06032E-02 0.00045 -2.40277E-05 0.00758 -1.84423E-04 0.01292  8.13228E-03 0.00892 ];
INF_S2                    (idx, [1:   8]) = [  1.11830E-02 0.00099 -2.74584E-06 0.06722 -1.14293E-04 0.04316  4.17905E-04 0.16074 ];
INF_S3                    (idx, [1:   8]) = [  3.11351E-03 0.00464 -4.99788E-07 0.22210 -4.24080E-05 0.15148  4.69998E-06 1.00000 ];
INF_S4                    (idx, [1:   8]) = [  1.54811E-03 0.00377 -3.44927E-07 0.59557 -1.95623E-05 0.20909  3.12400E-06 1.00000 ];
INF_S5                    (idx, [1:   8]) = [  5.21867E-04 0.01273 -2.22233E-07 0.27844 -8.08329E-06 0.74703  1.25205E-05 1.00000 ];
INF_S6                    (idx, [1:   8]) = [  2.29619E-04 0.04350 -2.46403E-09 1.00000 -7.66457E-06 0.62868  8.00107E-07 1.00000 ];
INF_S7                    (idx, [1:   8]) = [  7.87817E-05 0.07702 -6.76188E-08 1.00000 -4.46179E-06 0.56352  6.54502E-05 0.37797 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.43450E-01 9.5E-05  1.03696E-04 0.00344  2.25379E-03 0.00735  2.82506E-01 0.00016 ];
INF_SP1                   (idx, [1:   8]) = [  3.06034E-02 0.00045 -2.40277E-05 0.00758 -1.84423E-04 0.01292  8.13228E-03 0.00892 ];
INF_SP2                   (idx, [1:   8]) = [  1.11831E-02 0.00099 -2.74584E-06 0.06722 -1.14293E-04 0.04316  4.17905E-04 0.16074 ];
INF_SP3                   (idx, [1:   8]) = [  3.11351E-03 0.00464 -4.99788E-07 0.22210 -4.24080E-05 0.15148  4.69998E-06 1.00000 ];
INF_SP4                   (idx, [1:   8]) = [  1.54807E-03 0.00378 -3.44927E-07 0.59557 -1.95623E-05 0.20909  3.12400E-06 1.00000 ];
INF_SP5                   (idx, [1:   8]) = [  5.21905E-04 0.01272 -2.22233E-07 0.27844 -8.08329E-06 0.74703  1.25205E-05 1.00000 ];
INF_SP6                   (idx, [1:   8]) = [  2.29653E-04 0.04349 -2.46403E-09 1.00000 -7.66457E-06 0.62868  8.00107E-07 1.00000 ];
INF_SP7                   (idx, [1:   8]) = [  7.87874E-05 0.07725 -6.76188E-08 1.00000 -4.46179E-06 0.56352  6.54502E-05 0.37797 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  7.08606E-02 0.00211  4.12108E-02 0.00175 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  7.46726E-02 0.00223  3.50482E-02 0.00230 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  7.47044E-02 0.00198  3.50950E-02 0.00180 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  6.42727E-02 0.00215  6.34117E-02 0.00130 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  4.70415E+00 0.00211  8.08859E+00 0.00175 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  4.46402E+00 0.00223  9.51092E+00 0.00230 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  4.46210E+00 0.00198  9.49815E+00 0.00180 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  5.18633E+00 0.00215  5.25669E+00 0.00130 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  6.98374E-03 0.00563  2.27965E-04 0.03304  1.24418E-03 0.01451  1.17252E-03 0.01292  2.70921E-03 0.00795  1.16215E-03 0.01407  4.67719E-04 0.02044 ];
LAMBDA                    (idx, [1:  14]) = [  4.79860E-01 0.00815  1.33442E-02 0.00011  3.26738E-02 0.00011  1.20932E-01 6.0E-05  3.04356E-01 0.00015  8.55814E-01 0.00028  2.87498E+00 0.00043 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50812E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19145 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95768E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '5' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  4.92538E+05 0.00343  2.48993E+06 0.00229  5.80008E+06 0.00166  8.81672E+06 0.00217  1.01226E+07 0.00185  1.10100E+07 0.00170  6.29054E+06 0.00213  5.48323E+06 0.00204  9.59012E+06 0.00218  8.94873E+06 0.00212  1.00174E+07 0.00197  9.66602E+06 0.00162  9.69037E+06 0.00166  7.13836E+06 0.00151  5.42639E+06 0.00196  2.76594E+06 0.00161  7.02910E+05 0.00154  3.00154E+06 0.00172  2.86116E+06 0.00178  3.31977E+06 0.00155  1.58025E+06 0.00134  6.96742E+05 0.00198  3.41175E+05 0.00448  3.03000E+05 0.00472  2.61530E+05 0.00325  2.49341E+05 0.00400  3.93326E+05 0.00438  1.20927E+05 0.00327  1.57961E+05 0.00165  1.66254E+05 0.00568  8.98983E+04 0.00344  1.71984E+05 0.00307  1.14234E+05 0.00430  7.65511E+04 0.00434  1.13295E+04 0.00309  1.08362E+04 0.00698  1.12538E+04 0.00794  1.15495E+04 0.01135  1.15102E+04 0.00712  1.16130E+04 0.00719  1.18362E+04 0.01041  1.13338E+04 0.00518  2.15585E+04 0.00486  3.41680E+04 0.00591  4.34379E+04 0.00813  1.18609E+05 0.00270  1.33593E+05 0.00271  1.59905E+05 0.00062  1.05264E+05 0.00154  7.29268E+04 0.00179  5.37970E+04 0.00393  5.97785E+04 0.00291  1.06504E+05 0.00279  1.34430E+05 0.00208  2.28810E+05 0.00184  2.77047E+05 0.00123  2.90919E+05 0.00160  1.33152E+05 0.00230  7.62866E+04 0.00240  4.66176E+04 0.00358  3.61980E+04 0.00544  3.07459E+04 0.00363  2.20927E+04 0.00214  1.30244E+04 0.00204  1.04686E+04 0.00349  8.08918E+03 0.00413  5.78545E+03 0.00390  3.56087E+03 0.00202  1.74591E+03 0.00417  3.67252E+02 0.01504 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.56301E+00 0.00045 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.83857E+17 0.00183  7.56287E+15 0.00075 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  3.50796E-01 3.4E-05  5.73819E-01 0.00017 ];
INF_CAPT                  (idx, [1:   4]) = [  3.38056E-03 0.00039  5.70245E-02 0.00028 ];
INF_ABS                   (idx, [1:   4]) = [  7.21483E-03 0.00031  2.89197E-01 0.00033 ];
INF_FISS                  (idx, [1:   4]) = [  3.83426E-03 0.00030  2.32172E-01 0.00035 ];
INF_NSF                   (idx, [1:   4]) = [  9.47969E-03 0.00029  5.65734E-01 0.00035 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.47236E+00 7.2E-06  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02524E+02 6.1E-07  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  1.27662E-08 0.00048  1.82433E-06 0.00033 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.43586E-01 3.4E-05  2.84583E-01 0.00043 ];
INF_SCATT1                (idx, [1:   4]) = [  3.05859E-02 0.00049  7.78933E-03 0.01131 ];
INF_SCATT2                (idx, [1:   4]) = [  1.11687E-02 0.00055  2.47311E-04 0.31425 ];
INF_SCATT3                (idx, [1:   4]) = [  3.12998E-03 0.00228  6.88965E-05 1.00000 ];
INF_SCATT4                (idx, [1:   4]) = [  1.55158E-03 0.00408 -2.54836E-05 1.00000 ];
INF_SCATT5                (idx, [1:   4]) = [  5.19495E-04 0.00662 -6.61248E-05 0.39040 ];
INF_SCATT6                (idx, [1:   4]) = [  2.34096E-04 0.04747 -6.91523E-05 0.77287 ];
INF_SCATT7                (idx, [1:   4]) = [  8.71720E-05 0.05225  2.04105E-05 1.00000 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.43595E-01 3.4E-05  2.84583E-01 0.00043 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.05860E-02 0.00048  7.78933E-03 0.01131 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.11687E-02 0.00054  2.47311E-04 0.31425 ];
INF_SCATTP3               (idx, [1:   4]) = [  3.13005E-03 0.00227  6.88965E-05 1.00000 ];
INF_SCATTP4               (idx, [1:   4]) = [  1.55152E-03 0.00409 -2.54836E-05 1.00000 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.19500E-04 0.00661 -6.61248E-05 0.39040 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.34086E-04 0.04747 -6.91523E-05 0.77287 ];
INF_SCATTP7               (idx, [1:   4]) = [  8.72005E-05 0.05219  2.04105E-05 1.00000 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.83473E-01 5.6E-05  5.39554E-01 0.00028 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.17589E+00 5.6E-05  6.17794E-01 0.00028 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  7.20578E-03 0.00033  2.89197E-01 0.00033 ];
INF_REMXS                 (idx, [1:   4]) = [  7.31395E-03 0.00046  2.91477E-01 0.00075 ];

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

INF_S0                    (idx, [1:   8]) = [  3.43482E-01 3.3E-05  1.03284E-04 0.00361  2.24100E-03 0.00679  2.82342E-01 0.00045 ];
INF_S1                    (idx, [1:   8]) = [  3.06101E-02 0.00049 -2.41436E-05 0.00503 -1.82256E-04 0.02769  7.97159E-03 0.01141 ];
INF_S2                    (idx, [1:   8]) = [  1.11713E-02 0.00055 -2.64931E-06 0.03327 -1.00663E-04 0.07353  3.47974E-04 0.23952 ];
INF_S3                    (idx, [1:   8]) = [  3.13054E-03 0.00227 -5.60547E-07 0.25168 -4.42787E-05 0.10176  1.13175E-04 0.84476 ];
INF_S4                    (idx, [1:   8]) = [  1.55185E-03 0.00415 -2.69607E-07 0.50314 -1.78994E-05 0.22465 -7.58421E-06 1.00000 ];
INF_S5                    (idx, [1:   8]) = [  5.19714E-04 0.00690 -2.19530E-07 0.80911 -1.54740E-05 0.14774 -5.06508E-05 0.48637 ];
INF_S6                    (idx, [1:   8]) = [  2.34372E-04 0.04776 -2.75637E-07 0.44118 -8.05277E-06 0.28682 -6.10996E-05 0.88287 ];
INF_S7                    (idx, [1:   8]) = [  8.72117E-05 0.05249 -3.96795E-08 1.00000 -2.83458E-06 1.00000  2.32451E-05 1.00000 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.43491E-01 3.4E-05  1.03284E-04 0.00361  2.24100E-03 0.00679  2.82342E-01 0.00045 ];
INF_SP1                   (idx, [1:   8]) = [  3.06101E-02 0.00049 -2.41436E-05 0.00503 -1.82256E-04 0.02769  7.97159E-03 0.01141 ];
INF_SP2                   (idx, [1:   8]) = [  1.11714E-02 0.00055 -2.64931E-06 0.03327 -1.00663E-04 0.07353  3.47974E-04 0.23952 ];
INF_SP3                   (idx, [1:   8]) = [  3.13061E-03 0.00226 -5.60547E-07 0.25168 -4.42787E-05 0.10176  1.13175E-04 0.84476 ];
INF_SP4                   (idx, [1:   8]) = [  1.55179E-03 0.00416 -2.69607E-07 0.50314 -1.78994E-05 0.22465 -7.58421E-06 1.00000 ];
INF_SP5                   (idx, [1:   8]) = [  5.19720E-04 0.00688 -2.19530E-07 0.80911 -1.54740E-05 0.14774 -5.06508E-05 0.48637 ];
INF_SP6                   (idx, [1:   8]) = [  2.34362E-04 0.04776 -2.75637E-07 0.44118 -8.05277E-06 0.28682 -6.10996E-05 0.88287 ];
INF_SP7                   (idx, [1:   8]) = [  8.72402E-05 0.05243 -3.96795E-08 1.00000 -2.83458E-06 1.00000  2.32451E-05 1.00000 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  5.96933E-02 0.00192  3.23023E-02 0.00071 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  6.20558E-02 0.00174  2.40201E-02 0.00069 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  6.35856E-02 0.00196  3.06140E-02 0.00080 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  5.43019E-02 0.00205  5.38336E-02 0.00171 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  5.58419E+00 0.00192  1.03192E+01 0.00071 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  5.37157E+00 0.00174  1.38773E+01 0.00069 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  5.24236E+00 0.00196  1.08883E+01 0.00080 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  6.13862E+00 0.00205  6.19200E+00 0.00171 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  7.05714E-03 0.00534  2.22818E-04 0.02896  1.25755E-03 0.01258  1.18369E-03 0.01517  2.74099E-03 0.00800  1.16005E-03 0.01331  4.92039E-04 0.02072 ];
LAMBDA                    (idx, [1:  14]) = [  4.85996E-01 0.00808  1.33441E-02 9.1E-05  3.26724E-02 0.00011  1.20912E-01 6.6E-05  3.04370E-01 0.00016  8.55809E-01 0.00028  2.87643E+00 0.00041 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50812E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19145 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95768E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '6' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  4.94173E+05 0.00272  2.49745E+06 0.00141  5.80595E+06 0.00144  8.84359E+06 0.00154  1.01618E+07 0.00130  1.10372E+07 0.00145  6.30578E+06 0.00161  5.49811E+06 0.00164  9.60852E+06 0.00139  8.96817E+06 0.00105  1.00444E+07 0.00144  9.69081E+06 0.00158  9.71195E+06 0.00155  7.15970E+06 0.00164  5.44697E+06 0.00144  2.77904E+06 0.00161  7.06035E+05 0.00211  3.01403E+06 0.00164  2.87443E+06 0.00213  3.32704E+06 0.00210  1.58051E+06 0.00214  6.99006E+05 0.00391  3.41682E+05 0.00413  3.02940E+05 0.00191  2.62202E+05 0.00324  2.49819E+05 0.00350  3.93830E+05 0.00316  1.21006E+05 0.00402  1.57943E+05 0.00388  1.66744E+05 0.00517  9.07173E+04 0.00714  1.73512E+05 0.00277  1.15131E+05 0.00608  7.65609E+04 0.00350  1.12704E+04 0.00613  1.08892E+04 0.00863  1.12402E+04 0.00304  1.16394E+04 0.00696  1.17601E+04 0.00564  1.16840E+04 0.00282  1.21432E+04 0.01148  1.14331E+04 0.00699  2.15470E+04 0.00484  3.45909E+04 0.00465  4.37863E+04 0.00337  1.18250E+05 0.00203  1.35436E+05 0.00454  1.60183E+05 0.00355  1.05645E+05 0.00322  7.32795E+04 0.00231  5.41895E+04 0.00164  6.02162E+04 0.00320  1.07060E+05 0.00160  1.35154E+05 0.00245  2.30354E+05 0.00163  2.78446E+05 0.00093  2.92583E+05 0.00208  1.33664E+05 0.00244  7.65017E+04 0.00228  4.68060E+04 0.00241  3.61126E+04 0.00247  3.07971E+04 0.00334  2.19661E+04 0.00409  1.27652E+04 0.00698  1.03996E+04 0.00731  8.04158E+03 0.00654  5.72292E+03 0.00443  3.60756E+03 0.00850  1.71048E+03 0.00342  3.76500E+02 0.01265 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.56362E+00 0.00026 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.85175E+17 0.00157  7.60015E+15 0.00134 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  3.50816E-01 3.2E-05  5.73464E-01 0.00024 ];
INF_CAPT                  (idx, [1:   4]) = [  3.38288E-03 0.00047  5.69660E-02 0.00044 ];
INF_ABS                   (idx, [1:   4]) = [  7.21724E-03 0.00032  2.88852E-01 0.00046 ];
INF_FISS                  (idx, [1:   4]) = [  3.83436E-03 0.00021  2.31886E-01 0.00046 ];
INF_NSF                   (idx, [1:   4]) = [  9.47991E-03 0.00021  5.65036E-01 0.00046 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.47236E+00 6.5E-06  2.43670E+00 5.9E-09 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02524E+02 6.2E-07  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  1.27713E-08 0.00060  1.82249E-06 0.00037 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.43606E-01 3.0E-05  2.84777E-01 0.00040 ];
INF_SCATT1                (idx, [1:   4]) = [  3.05271E-02 0.00020  7.92113E-03 0.00806 ];
INF_SCATT2                (idx, [1:   4]) = [  1.11629E-02 0.00131  1.90471E-04 0.38958 ];
INF_SCATT3                (idx, [1:   4]) = [  3.11623E-03 0.00183 -7.34686E-05 0.75805 ];
INF_SCATT4                (idx, [1:   4]) = [  1.54176E-03 0.00250  8.23290E-05 0.57414 ];
INF_SCATT5                (idx, [1:   4]) = [  5.23723E-04 0.01028 -1.78092E-05 1.00000 ];
INF_SCATT6                (idx, [1:   4]) = [  2.30142E-04 0.02334 -7.79869E-05 0.35647 ];
INF_SCATT7                (idx, [1:   4]) = [  7.45816E-05 0.12789  1.29078E-05 1.00000 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.43615E-01 2.9E-05  2.84777E-01 0.00040 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.05272E-02 0.00020  7.92113E-03 0.00806 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.11629E-02 0.00131  1.90471E-04 0.38958 ];
INF_SCATTP3               (idx, [1:   4]) = [  3.11630E-03 0.00182 -7.34686E-05 0.75805 ];
INF_SCATTP4               (idx, [1:   4]) = [  1.54170E-03 0.00251  8.23290E-05 0.57414 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.23658E-04 0.01029 -1.78092E-05 1.00000 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.30079E-04 0.02338 -7.79869E-05 0.35647 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.45793E-05 0.12771  1.29078E-05 1.00000 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.83533E-01 6.7E-05  5.39184E-01 0.00036 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.17564E+00 6.7E-05  6.18218E-01 0.00036 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  7.20829E-03 0.00034  2.88852E-01 0.00046 ];
INF_REMXS                 (idx, [1:   4]) = [  7.31397E-03 0.00038  2.90937E-01 0.00044 ];

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

INF_S0                    (idx, [1:   8]) = [  3.43502E-01 3.0E-05  1.03782E-04 0.00402  2.24978E-03 0.01203  2.82527E-01 0.00049 ];
INF_S1                    (idx, [1:   8]) = [  3.05517E-02 0.00020 -2.45591E-05 0.01073 -1.93512E-04 0.04304  8.11464E-03 0.00858 ];
INF_S2                    (idx, [1:   8]) = [  1.11653E-02 0.00131 -2.41982E-06 0.07284 -1.10909E-04 0.05253  3.01380E-04 0.24520 ];
INF_S3                    (idx, [1:   8]) = [  3.11705E-03 0.00185 -8.18916E-07 0.26706 -4.39475E-05 0.11491 -2.95210E-05 1.00000 ];
INF_S4                    (idx, [1:   8]) = [  1.54197E-03 0.00250 -2.08096E-07 0.74325 -2.19525E-05 0.11802  1.04281E-04 0.46305 ];
INF_S5                    (idx, [1:   8]) = [  5.23869E-04 0.01032 -1.45504E-07 0.33701 -2.94439E-06 1.00000 -1.48648E-05 1.00000 ];
INF_S6                    (idx, [1:   8]) = [  2.30351E-04 0.02356 -2.08394E-07 0.54693 -1.37207E-06 1.00000 -7.66149E-05 0.39137 ];
INF_S7                    (idx, [1:   8]) = [  7.45815E-05 0.12829  4.12819E-11 1.00000 -1.86826E-06 1.00000  1.47761E-05 1.00000 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.43511E-01 2.9E-05  1.03782E-04 0.00402  2.24978E-03 0.01203  2.82527E-01 0.00049 ];
INF_SP1                   (idx, [1:   8]) = [  3.05518E-02 0.00020 -2.45591E-05 0.01073 -1.93512E-04 0.04304  8.11464E-03 0.00858 ];
INF_SP2                   (idx, [1:   8]) = [  1.11653E-02 0.00131 -2.41982E-06 0.07284 -1.10909E-04 0.05253  3.01380E-04 0.24520 ];
INF_SP3                   (idx, [1:   8]) = [  3.11712E-03 0.00185 -8.18916E-07 0.26706 -4.39475E-05 0.11491 -2.95210E-05 1.00000 ];
INF_SP4                   (idx, [1:   8]) = [  1.54191E-03 0.00251 -2.08096E-07 0.74325 -2.19525E-05 0.11802  1.04281E-04 0.46305 ];
INF_SP5                   (idx, [1:   8]) = [  5.23804E-04 0.01032 -1.45504E-07 0.33701 -2.94439E-06 1.00000 -1.48648E-05 1.00000 ];
INF_SP6                   (idx, [1:   8]) = [  2.30287E-04 0.02360 -2.08394E-07 0.54693 -1.37207E-06 1.00000 -7.66149E-05 0.39137 ];
INF_SP7                   (idx, [1:   8]) = [  7.45792E-05 0.12811  4.12819E-11 1.00000 -1.86826E-06 1.00000  1.47761E-05 1.00000 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  5.16478E-02 0.00155  2.66584E-02 0.00149 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  5.48208E-02 0.00154  2.43155E-02 0.00165 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  5.37434E-02 0.00156  1.99470E-02 0.00132 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  4.70863E-02 0.00158  4.70015E-02 0.00208 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  6.45404E+00 0.00155  1.25040E+01 0.00148 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  6.08048E+00 0.00154  1.37088E+01 0.00165 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  6.20237E+00 0.00156  1.67110E+01 0.00131 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  7.07927E+00 0.00158  7.09210E+00 0.00208 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  7.00474E-03 0.00525  2.33883E-04 0.02937  1.23604E-03 0.01300  1.16215E-03 0.01311  2.69676E-03 0.00859  1.17351E-03 0.01412  5.02402E-04 0.02174 ];
LAMBDA                    (idx, [1:  14]) = [  4.93070E-01 0.00885  1.33427E-02 8.8E-05  3.26719E-02 0.00012  1.20928E-01 6.6E-05  3.04299E-01 0.00015  8.55581E-01 0.00024  2.87327E+00 0.00044 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50812E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19144 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95768E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '7' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  4.92088E+05 0.00207  2.48871E+06 0.00180  5.80830E+06 0.00091  8.81919E+06 0.00099  1.01416E+07 0.00115  1.10201E+07 0.00108  6.29931E+06 0.00094  5.49608E+06 0.00076  9.60601E+06 0.00083  8.96507E+06 0.00077  1.00396E+07 0.00087  9.67871E+06 0.00117  9.69943E+06 0.00136  7.14498E+06 0.00140  5.44206E+06 0.00114  2.77366E+06 0.00100  7.04733E+05 0.00100  3.01029E+06 0.00092  2.87148E+06 0.00104  3.32337E+06 0.00040  1.57605E+06 0.00217  6.95749E+05 0.00183  3.41859E+05 0.00195  3.02966E+05 0.00357  2.62136E+05 0.00327  2.49210E+05 0.00387  3.93115E+05 0.00285  1.20952E+05 0.00369  1.58026E+05 0.00532  1.67088E+05 0.00322  9.03303E+04 0.00126  1.72404E+05 0.00296  1.14197E+05 0.00544  7.58725E+04 0.00502  1.12554E+04 0.00693  1.09435E+04 0.00662  1.10744E+04 0.00787  1.15447E+04 0.00572  1.15799E+04 0.00605  1.15395E+04 0.00378  1.20476E+04 0.00457  1.13139E+04 0.00714  2.15868E+04 0.00497  3.41415E+04 0.00331  4.37419E+04 0.00301  1.17954E+05 0.00217  1.35286E+05 0.00370  1.59756E+05 0.00370  1.06269E+05 0.00279  7.32382E+04 0.00247  5.39139E+04 0.00372  5.99256E+04 0.00326  1.06712E+05 0.00147  1.34350E+05 0.00168  2.29842E+05 0.00125  2.77444E+05 0.00165  2.91003E+05 0.00121  1.33561E+05 0.00184  7.62388E+04 0.00143  4.63825E+04 0.00140  3.62664E+04 0.00391  3.04888E+04 0.00298  2.19919E+04 0.00335  1.28562E+04 0.00243  1.04332E+04 0.00374  8.04946E+03 0.00340  5.78277E+03 0.00531  3.60091E+03 0.00414  1.73614E+03 0.00934  3.92855E+02 0.01444 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.56233E+00 0.00047 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.84560E+17 0.00100  7.57946E+15 0.00097 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  3.50823E-01 3.9E-05  5.73500E-01 0.00016 ];
INF_CAPT                  (idx, [1:   4]) = [  3.38021E-03 0.00015  5.69692E-02 0.00030 ];
INF_ABS                   (idx, [1:   4]) = [  7.21328E-03 0.00018  2.88885E-01 0.00030 ];
INF_FISS                  (idx, [1:   4]) = [  3.83307E-03 0.00023  2.31916E-01 0.00030 ];
INF_NSF                   (idx, [1:   4]) = [  9.47665E-03 0.00022  5.65110E-01 0.00030 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.47234E+00 8.6E-06  2.43670E+00 5.9E-09 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02524E+02 5.5E-07  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  1.27558E-08 0.00042  1.82262E-06 0.00029 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.43612E-01 4.0E-05  2.84397E-01 0.00034 ];
INF_SCATT1                (idx, [1:   4]) = [  3.05754E-02 0.00015  7.90229E-03 0.00895 ];
INF_SCATT2                (idx, [1:   4]) = [  1.11901E-02 0.00041  2.85249E-04 0.23986 ];
INF_SCATT3                (idx, [1:   4]) = [  3.11794E-03 0.00160  9.00353E-05 0.30227 ];
INF_SCATT4                (idx, [1:   4]) = [  1.55605E-03 0.00467  2.90621E-05 1.00000 ];
INF_SCATT5                (idx, [1:   4]) = [  5.25063E-04 0.00445  5.74583E-05 0.84492 ];
INF_SCATT6                (idx, [1:   4]) = [  2.36898E-04 0.01551 -5.73196E-06 1.00000 ];
INF_SCATT7                (idx, [1:   4]) = [  8.44274E-05 0.13285 -6.19749E-05 0.80174 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.43621E-01 4.0E-05  2.84397E-01 0.00034 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.05757E-02 0.00015  7.90229E-03 0.00895 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.11901E-02 0.00041  2.85249E-04 0.23986 ];
INF_SCATTP3               (idx, [1:   4]) = [  3.11791E-03 0.00159  9.00353E-05 0.30227 ];
INF_SCATTP4               (idx, [1:   4]) = [  1.55612E-03 0.00466  2.90621E-05 1.00000 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.25108E-04 0.00449  5.74583E-05 0.84492 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.36856E-04 0.01544 -5.73196E-06 1.00000 ];
INF_SCATTP7               (idx, [1:   4]) = [  8.44446E-05 0.13307 -6.19749E-05 0.80174 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.83524E-01 5.6E-05  5.39103E-01 0.00022 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.17568E+00 5.6E-05  6.18311E-01 0.00022 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  7.20423E-03 0.00019  2.88885E-01 0.00030 ];
INF_REMXS                 (idx, [1:   4]) = [  7.31371E-03 0.00051  2.91373E-01 0.00055 ];

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

INF_S0                    (idx, [1:   8]) = [  3.43509E-01 4.0E-05  1.03219E-04 0.00575  2.27018E-03 0.00608  2.82127E-01 0.00037 ];
INF_S1                    (idx, [1:   8]) = [  3.05998E-02 0.00015 -2.43588E-05 0.00609 -1.97490E-04 0.02584  8.09978E-03 0.00855 ];
INF_S2                    (idx, [1:   8]) = [  1.11925E-02 0.00040 -2.42722E-06 0.10929 -1.10592E-04 0.04799  3.95841E-04 0.17217 ];
INF_S3                    (idx, [1:   8]) = [  3.11866E-03 0.00157 -7.20112E-07 0.25082 -4.71058E-05 0.11725  1.37141E-04 0.21297 ];
INF_S4                    (idx, [1:   8]) = [  1.55628E-03 0.00464 -2.29991E-07 0.72527 -2.25688E-05 0.17365  5.16308E-05 1.00000 ];
INF_S5                    (idx, [1:   8]) = [  5.25102E-04 0.00439 -3.90442E-08 1.00000 -5.18145E-06 0.53038  6.26398E-05 0.77603 ];
INF_S6                    (idx, [1:   8]) = [  2.37155E-04 0.01569 -2.57004E-07 0.38996  3.08136E-06 1.00000 -8.81332E-06 1.00000 ];
INF_S7                    (idx, [1:   8]) = [  8.44001E-05 0.13266  2.72831E-08 1.00000 -5.33372E-06 0.72552 -5.66411E-05 0.86600 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.43518E-01 3.9E-05  1.03219E-04 0.00575  2.27018E-03 0.00608  2.82127E-01 0.00037 ];
INF_SP1                   (idx, [1:   8]) = [  3.06000E-02 0.00015 -2.43588E-05 0.00609 -1.97490E-04 0.02584  8.09978E-03 0.00855 ];
INF_SP2                   (idx, [1:   8]) = [  1.11926E-02 0.00040 -2.42722E-06 0.10929 -1.10592E-04 0.04799  3.95841E-04 0.17217 ];
INF_SP3                   (idx, [1:   8]) = [  3.11863E-03 0.00157 -7.20112E-07 0.25082 -4.71058E-05 0.11725  1.37141E-04 0.21297 ];
INF_SP4                   (idx, [1:   8]) = [  1.55635E-03 0.00464 -2.29991E-07 0.72527 -2.25688E-05 0.17365  5.16308E-05 1.00000 ];
INF_SP5                   (idx, [1:   8]) = [  5.25147E-04 0.00442 -3.90442E-08 1.00000 -5.18145E-06 0.53038  6.26398E-05 0.77603 ];
INF_SP6                   (idx, [1:   8]) = [  2.37113E-04 0.01562 -2.57004E-07 0.38996  3.08136E-06 1.00000 -8.81332E-06 1.00000 ];
INF_SP7                   (idx, [1:   8]) = [  8.44173E-05 0.13287  2.72831E-08 1.00000 -5.33372E-06 0.72552 -5.66411E-05 0.86600 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  4.53699E-02 0.00094  2.25639E-02 0.00103 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  4.76202E-02 0.00103  1.83725E-02 0.00108 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  4.76477E-02 0.00109  1.83770E-02 0.00131 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  4.14315E-02 0.00078  4.14751E-02 0.00127 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  7.34703E+00 0.00094  1.47730E+01 0.00103 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  6.99986E+00 0.00102  1.81431E+01 0.00108 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  6.99582E+00 0.00109  1.81387E+01 0.00131 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  8.04542E+00 0.00078  8.03701E+00 0.00127 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  6.99362E-03 0.00569  2.29453E-04 0.03158  1.23690E-03 0.01345  1.19252E-03 0.01494  2.69091E-03 0.00883  1.15677E-03 0.01388  4.87063E-04 0.01885 ];
LAMBDA                    (idx, [1:  14]) = [  4.86014E-01 0.00753  1.33428E-02 8.5E-05  3.26770E-02 0.00011  1.20919E-01 6.0E-05  3.04377E-01 0.00018  8.55622E-01 0.00028  2.87528E+00 0.00045 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50812E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19144 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95768E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '8' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  6.18030E+05 0.00152  2.94728E+06 0.00111  7.03562E+06 0.00061  9.56085E+06 0.00054  8.91667E+06 0.00025  9.65024E+06 0.00035  6.61406E+06 0.00061  5.79498E+06 0.00085  5.63442E+06 0.00037  4.72536E+06 0.00052  4.48757E+06 0.00054  4.27061E+06 0.00064  4.12606E+06 0.00029  3.90950E+06 0.00058  3.76051E+06 0.00044  2.82805E+06 0.00083  1.74100E+06 0.00056  4.29608E+06 0.00048  3.35216E+06 0.00078  6.38724E+06 0.00033  6.18069E+06 0.00048  4.51147E+06 0.00046  2.94933E+06 0.00069  3.53484E+06 0.00071  3.46779E+06 0.00049  2.98531E+06 0.00052  5.51918E+06 0.00053  1.17272E+06 0.00053  1.45442E+06 0.00060  1.30919E+06 0.00036  7.62709E+05 0.00060  1.31467E+06 0.00040  8.93699E+05 0.00092  7.69432E+05 0.00039  1.50057E+05 0.00122  1.47907E+05 0.00195  1.52205E+05 0.00132  1.56692E+05 0.00082  1.54811E+05 0.00164  1.52958E+05 0.00122  1.57576E+05 0.00061  1.48869E+05 0.00102  2.82180E+05 0.00096  4.56822E+05 0.00062  6.01526E+05 0.00055  1.82655E+06 0.00079  2.80540E+06 0.00037  4.85398E+06 0.00048  4.37109E+06 0.00048  3.66347E+06 0.00043  3.01552E+06 0.00052  3.57128E+06 0.00047  6.60623E+06 0.00035  8.41571E+06 0.00028  1.46613E+07 0.00036  1.91986E+07 0.00037  2.36401E+07 0.00034  1.29349E+07 0.00049  8.46338E+06 0.00039  5.72628E+06 0.00078  4.89810E+06 0.00046  4.61842E+06 0.00052  3.71950E+06 0.00046  2.42061E+06 0.00048  2.18339E+06 0.00044  1.90308E+06 0.00050  1.57855E+06 0.00070  1.20123E+06 0.00039  7.67744E+05 0.00050  2.70780E+05 0.00082 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  5.33509E+17 0.00049  5.47118E+17 0.00047 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.54682E-01 5.9E-05  8.93380E-01 1.1E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  4.94718E-04 0.00019  1.04628E-02 2.4E-05 ];
INF_ABS                   (idx, [1:   4]) = [  4.94718E-04 0.00019  1.04628E-02 2.4E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  8.45662E-08 0.00025  2.24365E-06 2.4E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.54187E-01 5.9E-05  8.82916E-01 1.1E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  2.39976E-01 7.9E-05  2.77632E-01 2.4E-05 ];
INF_SCATT2                (idx, [1:   4]) = [  9.22482E-02 0.00018  1.01185E-01 8.5E-05 ];
INF_SCATT3                (idx, [1:   4]) = [  4.94606E-03 0.00112  4.13863E-02 0.00040 ];
INF_SCATT4                (idx, [1:   4]) = [ -1.27441E-02 0.00071  1.99999E-02 0.00065 ];
INF_SCATT5                (idx, [1:   4]) = [ -1.92001E-03 0.00262  1.14229E-02 0.00091 ];
INF_SCATT6                (idx, [1:   4]) = [  4.00400E-03 0.00194  7.41429E-03 0.00116 ];
INF_SCATT7                (idx, [1:   4]) = [  7.49015E-04 0.01258  5.26752E-03 0.00257 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.54187E-01 5.9E-05  8.82916E-01 1.1E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.39976E-01 8.0E-05  2.77632E-01 2.4E-05 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.22482E-02 0.00018  1.01185E-01 8.5E-05 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.94606E-03 0.00112  4.13863E-02 0.00040 ];
INF_SCATTP4               (idx, [1:   4]) = [ -1.27441E-02 0.00071  1.99999E-02 0.00065 ];
INF_SCATTP5               (idx, [1:   4]) = [ -1.92001E-03 0.00262  1.14229E-02 0.00091 ];
INF_SCATTP6               (idx, [1:   4]) = [  4.00400E-03 0.00194  7.41429E-03 0.00116 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.49016E-04 0.01258  5.26752E-03 0.00257 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.30643E-01 0.00010  5.71470E-01 1.5E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.44524E+00 0.00010  5.83291E-01 1.5E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  4.94688E-04 0.00019  1.04628E-02 2.4E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  2.38979E-02 0.00017  1.26053E-02 0.00024 ];

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

INF_CHIT                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.30785E-01 6.4E-05  2.34027E-02 0.00018  2.14183E-03 0.00128  8.80775E-01 1.3E-05 ];
INF_S1                    (idx, [1:   8]) = [  2.32028E-01 8.4E-05  7.94851E-03 0.00040  9.24970E-04 0.00137  2.76707E-01 2.6E-05 ];
INF_S2                    (idx, [1:   8]) = [  9.31018E-02 0.00016 -8.53567E-04 0.00342  3.50350E-04 0.00253  1.00834E-01 9.1E-05 ];
INF_S3                    (idx, [1:   8]) = [  7.31388E-03 0.00074 -2.36782E-03 0.00097  5.78216E-05 0.01260  4.13285E-02 0.00041 ];
INF_S4                    (idx, [1:   8]) = [ -1.15581E-02 0.00080 -1.18595E-03 0.00152 -4.68213E-05 0.00690  2.00467E-02 0.00064 ];
INF_S5                    (idx, [1:   8]) = [ -1.58118E-03 0.00320 -3.38833E-04 0.00689 -6.43228E-05 0.00696  1.14872E-02 0.00090 ];
INF_S6                    (idx, [1:   8]) = [  4.11716E-03 0.00199 -1.13165E-04 0.01497 -5.34091E-05 0.00852  7.46770E-03 0.00110 ];
INF_S7                    (idx, [1:   8]) = [  8.23827E-04 0.01351 -7.48116E-05 0.02629 -3.99775E-05 0.00836  5.30750E-03 0.00254 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.30785E-01 6.4E-05  2.34027E-02 0.00018  2.14183E-03 0.00128  8.80775E-01 1.3E-05 ];
INF_SP1                   (idx, [1:   8]) = [  2.32028E-01 8.4E-05  7.94851E-03 0.00040  9.24970E-04 0.00137  2.76707E-01 2.6E-05 ];
INF_SP2                   (idx, [1:   8]) = [  9.31018E-02 0.00016 -8.53567E-04 0.00342  3.50350E-04 0.00253  1.00834E-01 9.1E-05 ];
INF_SP3                   (idx, [1:   8]) = [  7.31388E-03 0.00074 -2.36782E-03 0.00097  5.78216E-05 0.01260  4.13285E-02 0.00041 ];
INF_SP4                   (idx, [1:   8]) = [ -1.15581E-02 0.00080 -1.18595E-03 0.00152 -4.68213E-05 0.00690  2.00467E-02 0.00064 ];
INF_SP5                   (idx, [1:   8]) = [ -1.58118E-03 0.00320 -3.38833E-04 0.00689 -6.43228E-05 0.00696  1.14872E-02 0.00090 ];
INF_SP6                   (idx, [1:   8]) = [  4.11716E-03 0.00199 -1.13165E-04 0.01497 -5.34091E-05 0.00852  7.46770E-03 0.00110 ];
INF_SP7                   (idx, [1:   8]) = [  8.23827E-04 0.01352 -7.48116E-05 0.02629 -3.99775E-05 0.00836  5.30750E-03 0.00254 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  3.96011E-02 0.00050  7.64423E-01 0.00038 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  4.05419E-02 0.00042  1.02476E+00 0.00110 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  4.05512E-02 0.00058  1.02282E+00 0.00118 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  3.78365E-02 0.00059  5.07359E-01 0.00059 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  8.41729E+00 0.00049  4.36059E-01 0.00038 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  8.22195E+00 0.00042  3.25282E-01 0.00110 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  8.22006E+00 0.00058  3.25897E-01 0.00118 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  8.80985E+00 0.00059  6.56998E-01 0.00059 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50813E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19142 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95768E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '9' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  2.54932E+05 0.00557  1.21319E+06 0.00206  2.86936E+06 0.00096  3.86876E+06 0.00146  3.57002E+06 0.00059  3.84533E+06 0.00107  2.63244E+06 0.00124  2.30354E+06 0.00118  2.23510E+06 0.00093  1.87676E+06 0.00107  1.78018E+06 0.00092  1.69832E+06 0.00094  1.64150E+06 0.00109  1.55840E+06 0.00145  1.49889E+06 0.00101  1.12522E+06 0.00146  6.91454E+05 0.00110  1.70461E+06 0.00145  1.32848E+06 0.00137  2.53492E+06 0.00113  2.44730E+06 0.00133  1.78727E+06 0.00078  1.16677E+06 0.00093  1.39845E+06 0.00123  1.37044E+06 0.00117  1.17885E+06 0.00126  2.17488E+06 0.00083  4.61410E+05 0.00156  5.71707E+05 0.00141  5.14518E+05 0.00135  2.99676E+05 0.00114  5.15727E+05 0.00150  3.50430E+05 0.00110  3.01233E+05 0.00162  5.84261E+04 0.00290  5.81282E+04 0.00364  5.95643E+04 0.00189  6.12742E+04 0.00131  6.05160E+04 0.00241  5.98993E+04 0.00379  6.14020E+04 0.00349  5.80856E+04 0.00269  1.10191E+05 0.00342  1.78184E+05 0.00187  2.33596E+05 0.00141  7.04363E+05 0.00227  1.05876E+06 0.00099  1.79478E+06 0.00137  1.60462E+06 0.00131  1.33828E+06 0.00133  1.09851E+06 0.00141  1.30042E+06 0.00113  2.40173E+06 0.00100  3.05286E+06 0.00119  5.31460E+06 0.00083  6.94726E+06 0.00065  8.54663E+06 0.00092  4.67542E+06 0.00074  3.05747E+06 0.00091  2.06938E+06 0.00096  1.76894E+06 0.00087  1.66792E+06 0.00104  1.34334E+06 0.00111  8.74672E+05 0.00088  7.88507E+05 0.00075  6.87315E+05 0.00103  5.70014E+05 0.00117  4.34222E+05 0.00101  2.77089E+05 0.00136  9.79713E+04 0.00065 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  2.12372E+17 0.00103  1.98447E+17 0.00085 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.52828E-01 5.6E-05  8.92791E-01 1.7E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  4.90530E-04 0.00021  1.04446E-02 4.1E-05 ];
INF_ABS                   (idx, [1:   4]) = [  4.90530E-04 0.00021  1.04446E-02 4.1E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  8.33347E-08 0.00033  2.23974E-06 4.1E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.52339E-01 5.6E-05  8.82342E-01 2.3E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  2.39067E-01 7.6E-05  2.77756E-01 6.1E-05 ];
INF_SCATT2                (idx, [1:   4]) = [  9.19143E-02 0.00024  1.01242E-01 0.00028 ];
INF_SCATT3                (idx, [1:   4]) = [  4.89455E-03 0.00117  4.14129E-02 0.00070 ];
INF_SCATT4                (idx, [1:   4]) = [ -1.27094E-02 0.00130  1.99598E-02 0.00167 ];
INF_SCATT5                (idx, [1:   4]) = [ -1.90789E-03 0.00435  1.13605E-02 0.00234 ];
INF_SCATT6                (idx, [1:   4]) = [  3.97364E-03 0.00272  7.36577E-03 0.00285 ];
INF_SCATT7                (idx, [1:   4]) = [  7.12200E-04 0.01901  5.21837E-03 0.00234 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.52339E-01 5.6E-05  8.82342E-01 2.3E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.39067E-01 7.6E-05  2.77756E-01 6.1E-05 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.19143E-02 0.00024  1.01242E-01 0.00028 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.89456E-03 0.00117  4.14129E-02 0.00070 ];
INF_SCATTP4               (idx, [1:   4]) = [ -1.27094E-02 0.00130  1.99598E-02 0.00167 ];
INF_SCATTP5               (idx, [1:   4]) = [ -1.90789E-03 0.00435  1.13605E-02 0.00234 ];
INF_SCATTP6               (idx, [1:   4]) = [  3.97365E-03 0.00272  7.36577E-03 0.00285 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.12198E-04 0.01901  5.21837E-03 0.00234 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.29224E-01 0.00013  5.70786E-01 1.7E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.45418E+00 0.00013  5.83990E-01 1.7E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  4.90515E-04 0.00022  1.04446E-02 4.1E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  2.34462E-02 0.00016  1.26442E-02 0.00071 ];

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

INF_CHIT                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.29382E-01 6.3E-05  2.29564E-02 0.00018  2.19501E-03 0.00042  8.80147E-01 2.4E-05 ];
INF_S1                    (idx, [1:   8]) = [  2.31278E-01 8.1E-05  7.78892E-03 0.00083  9.50052E-04 0.00145  2.76806E-01 6.0E-05 ];
INF_S2                    (idx, [1:   8]) = [  9.27750E-02 0.00022 -8.60750E-04 0.00326  3.63262E-04 0.00299  1.00879E-01 0.00027 ];
INF_S3                    (idx, [1:   8]) = [  7.23669E-03 0.00064 -2.34214E-03 0.00139  6.15484E-05 0.01952  4.13514E-02 0.00068 ];
INF_S4                    (idx, [1:   8]) = [ -1.15462E-02 0.00133 -1.16323E-03 0.00133 -4.65635E-05 0.02561  2.00063E-02 0.00164 ];
INF_S5                    (idx, [1:   8]) = [ -1.58364E-03 0.00665 -3.24253E-04 0.00863 -6.55318E-05 0.01471  1.14260E-02 0.00236 ];
INF_S6                    (idx, [1:   8]) = [  4.08056E-03 0.00302 -1.06920E-04 0.01804 -5.56850E-05 0.01587  7.42146E-03 0.00288 ];
INF_S7                    (idx, [1:   8]) = [  7.84735E-04 0.01632 -7.25344E-05 0.02144 -4.05658E-05 0.00850  5.25894E-03 0.00237 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.29382E-01 6.3E-05  2.29564E-02 0.00018  2.19501E-03 0.00042  8.80147E-01 2.4E-05 ];
INF_SP1                   (idx, [1:   8]) = [  2.31278E-01 8.1E-05  7.78892E-03 0.00083  9.50052E-04 0.00145  2.76806E-01 6.0E-05 ];
INF_SP2                   (idx, [1:   8]) = [  9.27750E-02 0.00022 -8.60750E-04 0.00326  3.63262E-04 0.00299  1.00879E-01 0.00027 ];
INF_SP3                   (idx, [1:   8]) = [  7.23670E-03 0.00064 -2.34214E-03 0.00139  6.15484E-05 0.01952  4.13514E-02 0.00068 ];
INF_SP4                   (idx, [1:   8]) = [ -1.15462E-02 0.00133 -1.16323E-03 0.00133 -4.65635E-05 0.02561  2.00063E-02 0.00164 ];
INF_SP5                   (idx, [1:   8]) = [ -1.58364E-03 0.00665 -3.24253E-04 0.00863 -6.55318E-05 0.01471  1.14260E-02 0.00236 ];
INF_SP6                   (idx, [1:   8]) = [  4.08057E-03 0.00302 -1.06920E-04 0.01804 -5.56850E-05 0.01587  7.42146E-03 0.00288 ];
INF_SP7                   (idx, [1:   8]) = [  7.84732E-04 0.01632 -7.25344E-05 0.02144 -4.05658E-05 0.00850  5.25894E-03 0.00237 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  1.43511E-02 0.00091  1.77051E-01 0.00102 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  1.44081E-02 0.00101  1.79263E-01 0.00101 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  1.45714E-02 0.00098  2.32711E-01 0.00133 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  1.40824E-02 0.00079  1.41469E-01 0.00113 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  2.32271E+01 0.00091  1.88271E+00 0.00102 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  2.31352E+01 0.00101  1.85947E+00 0.00101 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  2.28760E+01 0.00098  1.43240E+00 0.00133 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  2.36703E+01 0.00079  2.35625E+00 0.00113 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50813E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19143 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95768E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  2])  = '10' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  2.55780E+05 0.00299  1.20940E+06 0.00205  2.86497E+06 0.00161  3.86966E+06 0.00219  3.56559E+06 0.00152  3.83664E+06 0.00129  2.62519E+06 0.00179  2.29718E+06 0.00208  2.23113E+06 0.00156  1.87095E+06 0.00196  1.77740E+06 0.00193  1.69313E+06 0.00143  1.63855E+06 0.00157  1.55404E+06 0.00161  1.49529E+06 0.00228  1.12219E+06 0.00174  6.89876E+05 0.00208  1.69992E+06 0.00205  1.32772E+06 0.00178  2.53128E+06 0.00183  2.44544E+06 0.00164  1.78564E+06 0.00170  1.16568E+06 0.00175  1.39509E+06 0.00123  1.36803E+06 0.00188  1.17815E+06 0.00165  2.17443E+06 0.00111  4.61344E+05 0.00282  5.71388E+05 0.00077  5.14708E+05 0.00270  2.99249E+05 0.00225  5.15379E+05 0.00191  3.50289E+05 0.00190  3.01934E+05 0.00242  5.86734E+04 0.00493  5.80548E+04 0.00243  5.94205E+04 0.00351  6.11501E+04 0.00250  6.03329E+04 0.00221  5.97595E+04 0.00152  6.14453E+04 0.00246  5.80085E+04 0.00238  1.09973E+05 0.00210  1.78058E+05 0.00149  2.33914E+05 0.00149  7.04828E+05 0.00146  1.05794E+06 0.00237  1.79428E+06 0.00133  1.60091E+06 0.00171  1.33609E+06 0.00174  1.09789E+06 0.00214  1.29856E+06 0.00161  2.39749E+06 0.00127  3.05085E+06 0.00158  5.30433E+06 0.00154  6.93887E+06 0.00173  8.53434E+06 0.00158  4.66906E+06 0.00173  3.05286E+06 0.00179  2.06630E+06 0.00153  1.76689E+06 0.00173  1.66545E+06 0.00159  1.34194E+06 0.00161  8.73928E+05 0.00159  7.87633E+05 0.00187  6.86793E+05 0.00187  5.69586E+05 0.00202  4.33508E+05 0.00208  2.77188E+05 0.00238  9.77415E+04 0.00208 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  2.12043E+17 0.00160  1.98183E+17 0.00159 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.52798E-01 8.5E-05  8.92810E-01 2.7E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  4.91003E-04 0.00031  1.04450E-02 6.5E-05 ];
INF_ABS                   (idx, [1:   4]) = [  4.91003E-04 0.00031  1.04450E-02 6.5E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  8.34295E-08 0.00029  2.23983E-06 6.5E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.52306E-01 8.6E-05  8.82363E-01 3.1E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  2.39068E-01 7.2E-05  2.77680E-01 6.4E-05 ];
INF_SCATT2                (idx, [1:   4]) = [  9.18898E-02 0.00017  1.01229E-01 0.00025 ];
INF_SCATT3                (idx, [1:   4]) = [  4.85133E-03 0.00235  4.14152E-02 0.00056 ];
INF_SCATT4                (idx, [1:   4]) = [ -1.27214E-02 0.00073  1.99514E-02 0.00068 ];
INF_SCATT5                (idx, [1:   4]) = [ -1.89702E-03 0.00815  1.13650E-02 0.00117 ];
INF_SCATT6                (idx, [1:   4]) = [  3.98666E-03 0.00273  7.36839E-03 0.00292 ];
INF_SCATT7                (idx, [1:   4]) = [  7.29661E-04 0.01820  5.23851E-03 0.00525 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.52306E-01 8.6E-05  8.82363E-01 3.1E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.39068E-01 7.2E-05  2.77680E-01 6.4E-05 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.18898E-02 0.00017  1.01229E-01 0.00025 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.85133E-03 0.00234  4.14152E-02 0.00056 ];
INF_SCATTP4               (idx, [1:   4]) = [ -1.27214E-02 0.00073  1.99514E-02 0.00068 ];
INF_SCATTP5               (idx, [1:   4]) = [ -1.89702E-03 0.00815  1.13650E-02 0.00117 ];
INF_SCATTP6               (idx, [1:   4]) = [  3.98666E-03 0.00273  7.36839E-03 0.00292 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.29662E-04 0.01820  5.23851E-03 0.00525 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.29214E-01 0.00018  5.70879E-01 3.0E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.45424E+00 0.00018  5.83895E-01 3.0E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  4.90975E-04 0.00030  1.04450E-02 6.5E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  2.34690E-02 0.00026  1.26431E-02 0.00074 ];

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

INF_CHIT                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.29329E-01 8.6E-05  2.29769E-02 0.00028  2.19611E-03 0.00175  8.80167E-01 3.4E-05 ];
INF_S1                    (idx, [1:   8]) = [  2.31272E-01 7.1E-05  7.79603E-03 0.00073  9.51906E-04 0.00241  2.76729E-01 6.9E-05 ];
INF_S2                    (idx, [1:   8]) = [  9.27513E-02 0.00020 -8.61524E-04 0.00330  3.66352E-04 0.00645  1.00862E-01 0.00024 ];
INF_S3                    (idx, [1:   8]) = [  7.19674E-03 0.00155 -2.34541E-03 0.00112  6.36555E-05 0.02558  4.13516E-02 0.00054 ];
INF_S4                    (idx, [1:   8]) = [ -1.15536E-02 0.00083 -1.16776E-03 0.00281 -4.41494E-05 0.02370  1.99956E-02 0.00068 ];
INF_S5                    (idx, [1:   8]) = [ -1.57372E-03 0.01076 -3.23301E-04 0.01087 -6.39028E-05 0.00838  1.14289E-02 0.00116 ];
INF_S6                    (idx, [1:   8]) = [  4.09070E-03 0.00270 -1.04038E-04 0.01148 -5.37304E-05 0.01516  7.42212E-03 0.00286 ];
INF_S7                    (idx, [1:   8]) = [  8.02563E-04 0.01435 -7.29016E-05 0.03192 -4.05648E-05 0.01989  5.27907E-03 0.00528 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.29329E-01 8.6E-05  2.29769E-02 0.00028  2.19611E-03 0.00175  8.80167E-01 3.4E-05 ];
INF_SP1                   (idx, [1:   8]) = [  2.31272E-01 7.1E-05  7.79603E-03 0.00073  9.51906E-04 0.00241  2.76729E-01 6.9E-05 ];
INF_SP2                   (idx, [1:   8]) = [  9.27513E-02 0.00020 -8.61524E-04 0.00330  3.66352E-04 0.00645  1.00862E-01 0.00024 ];
INF_SP3                   (idx, [1:   8]) = [  7.19674E-03 0.00155 -2.34541E-03 0.00112  6.36555E-05 0.02558  4.13516E-02 0.00054 ];
INF_SP4                   (idx, [1:   8]) = [ -1.15536E-02 0.00083 -1.16776E-03 0.00281 -4.41494E-05 0.02370  1.99956E-02 0.00068 ];
INF_SP5                   (idx, [1:   8]) = [ -1.57372E-03 0.01076 -3.23301E-04 0.01087 -6.39028E-05 0.00838  1.14289E-02 0.00116 ];
INF_SP6                   (idx, [1:   8]) = [  4.09070E-03 0.00270 -1.04038E-04 0.01148 -5.37304E-05 0.01516  7.42212E-03 0.00286 ];
INF_SP7                   (idx, [1:   8]) = [  8.02564E-04 0.01434 -7.29016E-05 0.03192 -4.05648E-05 0.01989  5.27907E-03 0.00528 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  1.31528E-02 0.00150  1.29859E-01 0.00145 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  1.31928E-02 0.00162  1.52784E-01 0.00160 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  1.30617E-02 0.00153  1.27542E-01 0.00133 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  1.32047E-02 0.00137  1.14729E-01 0.00169 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  2.53435E+01 0.00150  2.56691E+00 0.00145 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  2.52665E+01 0.00162  2.18176E+00 0.00160 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  2.55201E+01 0.00152  2.61354E+00 0.00133 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  2.52438E+01 0.00136  2.90542E+00 0.00168 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50813E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19141 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95767E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  2])  = '11' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  2.57069E+05 0.00334  1.21273E+06 0.00238  2.86969E+06 0.00191  3.87783E+06 0.00229  3.57686E+06 0.00226  3.84803E+06 0.00227  2.63283E+06 0.00226  2.30221E+06 0.00169  2.23454E+06 0.00243  1.87350E+06 0.00243  1.78158E+06 0.00220  1.69650E+06 0.00239  1.64070E+06 0.00184  1.55934E+06 0.00236  1.49912E+06 0.00177  1.12360E+06 0.00220  6.92004E+05 0.00229  1.70338E+06 0.00226  1.32964E+06 0.00196  2.53780E+06 0.00207  2.45050E+06 0.00215  1.78725E+06 0.00222  1.16828E+06 0.00221  1.39924E+06 0.00147  1.37021E+06 0.00161  1.17895E+06 0.00231  2.17594E+06 0.00163  4.62048E+05 0.00244  5.72235E+05 0.00224  5.14588E+05 0.00224  2.99768E+05 0.00174  5.16580E+05 0.00228  3.49621E+05 0.00207  3.01471E+05 0.00134  5.86911E+04 0.00315  5.78612E+04 0.00244  5.95021E+04 0.00220  6.10616E+04 0.00139  6.06427E+04 0.00413  5.97796E+04 0.00235  6.15409E+04 0.00334  5.82057E+04 0.00143  1.10516E+05 0.00180  1.78578E+05 0.00271  2.33886E+05 0.00275  7.05152E+05 0.00253  1.05854E+06 0.00197  1.79526E+06 0.00149  1.60247E+06 0.00155  1.33889E+06 0.00213  1.09910E+06 0.00233  1.29919E+06 0.00208  2.39749E+06 0.00166  3.05012E+06 0.00196  5.30759E+06 0.00187  6.94101E+06 0.00189  8.53912E+06 0.00181  4.66873E+06 0.00171  3.05343E+06 0.00143  2.06616E+06 0.00142  1.76748E+06 0.00128  1.66584E+06 0.00185  1.34209E+06 0.00151  8.73400E+05 0.00156  7.87341E+05 0.00149  6.87146E+05 0.00126  5.69566E+05 0.00151  4.33795E+05 0.00172  2.76801E+05 0.00181  9.77872E+04 0.00165 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  2.12475E+17 0.00195  1.98250E+17 0.00165 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.52767E-01 0.00012  8.92770E-01 3.6E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  4.90701E-04 0.00039  1.04439E-02 8.3E-05 ];
INF_ABS                   (idx, [1:   4]) = [  4.90701E-04 0.00039  1.04439E-02 8.3E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  8.33444E-08 0.00031  2.23960E-06 8.3E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.52277E-01 0.00012  8.82331E-01 4.0E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  2.39066E-01 0.00011  2.77764E-01 0.00011 ];
INF_SCATT2                (idx, [1:   4]) = [  9.18944E-02 0.00020  1.01227E-01 0.00032 ];
INF_SCATT3                (idx, [1:   4]) = [  4.88407E-03 0.00369  4.13782E-02 0.00056 ];
INF_SCATT4                (idx, [1:   4]) = [ -1.26996E-02 0.00103  1.99708E-02 0.00104 ];
INF_SCATT5                (idx, [1:   4]) = [ -1.90681E-03 0.00402  1.13786E-02 0.00229 ];
INF_SCATT6                (idx, [1:   4]) = [  4.01043E-03 0.00275  7.37251E-03 0.00151 ];
INF_SCATT7                (idx, [1:   4]) = [  7.59020E-04 0.01154  5.24727E-03 0.00410 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.52277E-01 0.00012  8.82331E-01 4.0E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.39066E-01 0.00011  2.77764E-01 0.00011 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.18943E-02 0.00020  1.01227E-01 0.00032 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.88408E-03 0.00369  4.13782E-02 0.00056 ];
INF_SCATTP4               (idx, [1:   4]) = [ -1.26996E-02 0.00103  1.99708E-02 0.00104 ];
INF_SCATTP5               (idx, [1:   4]) = [ -1.90681E-03 0.00402  1.13786E-02 0.00229 ];
INF_SCATTP6               (idx, [1:   4]) = [  4.01043E-03 0.00275  7.37251E-03 0.00151 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.59020E-04 0.01154  5.24727E-03 0.00410 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.29180E-01 0.00024  5.70750E-01 5.0E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.45446E+00 0.00024  5.84027E-01 5.0E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  4.90676E-04 0.00039  1.04439E-02 8.3E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  2.34414E-02 0.00024  1.26351E-02 0.00081 ];

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

INF_CHIT                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.29325E-01 0.00012  2.29517E-02 0.00023  2.19571E-03 0.00134  8.80135E-01 4.2E-05 ];
INF_S1                    (idx, [1:   8]) = [  2.31283E-01 0.00010  7.78296E-03 0.00053  9.51617E-04 0.00234  2.76812E-01 0.00011 ];
INF_S2                    (idx, [1:   8]) = [  9.27559E-02 0.00016 -8.61530E-04 0.00546  3.65208E-04 0.00399  1.00861E-01 0.00033 ];
INF_S3                    (idx, [1:   8]) = [  7.22122E-03 0.00203 -2.33715E-03 0.00164  6.29520E-05 0.02382  4.13152E-02 0.00060 ];
INF_S4                    (idx, [1:   8]) = [ -1.15390E-02 0.00095 -1.16058E-03 0.00400 -4.45429E-05 0.01176  2.00154E-02 0.00104 ];
INF_S5                    (idx, [1:   8]) = [ -1.57887E-03 0.00548 -3.27932E-04 0.00685 -6.45783E-05 0.02146  1.14432E-02 0.00229 ];
INF_S6                    (idx, [1:   8]) = [  4.12304E-03 0.00242 -1.12610E-04 0.00972 -5.45083E-05 0.02395  7.42702E-03 0.00154 ];
INF_S7                    (idx, [1:   8]) = [  8.31380E-04 0.00971 -7.23602E-05 0.02461 -4.01571E-05 0.01816  5.28743E-03 0.00413 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.29325E-01 0.00012  2.29517E-02 0.00023  2.19571E-03 0.00134  8.80135E-01 4.2E-05 ];
INF_SP1                   (idx, [1:   8]) = [  2.31283E-01 0.00010  7.78296E-03 0.00053  9.51617E-04 0.00234  2.76812E-01 0.00011 ];
INF_SP2                   (idx, [1:   8]) = [  9.27559E-02 0.00016 -8.61530E-04 0.00546  3.65208E-04 0.00399  1.00861E-01 0.00033 ];
INF_SP3                   (idx, [1:   8]) = [  7.22122E-03 0.00203 -2.33715E-03 0.00164  6.29520E-05 0.02382  4.13152E-02 0.00060 ];
INF_SP4                   (idx, [1:   8]) = [ -1.15390E-02 0.00095 -1.16058E-03 0.00400 -4.45429E-05 0.01176  2.00154E-02 0.00104 ];
INF_SP5                   (idx, [1:   8]) = [ -1.57887E-03 0.00548 -3.27932E-04 0.00685 -6.45783E-05 0.02146  1.14432E-02 0.00229 ];
INF_SP6                   (idx, [1:   8]) = [  4.12304E-03 0.00242 -1.12610E-04 0.00972 -5.45083E-05 0.02395  7.42702E-03 0.00154 ];
INF_SP7                   (idx, [1:   8]) = [  8.31380E-04 0.00971 -7.23602E-05 0.02461 -4.01571E-05 0.01816  5.28743E-03 0.00413 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  1.21778E-02 0.00187  1.02663E-01 0.00140 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  1.20366E-02 0.00189  1.06061E-01 0.00182 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  1.20371E-02 0.00182  1.05917E-01 0.00118 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  1.24698E-02 0.00191  9.65993E-02 0.00128 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  2.73727E+01 0.00187  3.24690E+00 0.00141 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  2.76937E+01 0.00189  3.14288E+00 0.00182 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  2.76926E+01 0.00182  3.14713E+00 0.00118 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  2.67317E+01 0.00191  3.45070E+00 0.00129 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50813E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19141 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95767E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  2])  = '12' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  2.57830E+05 0.00174  1.21167E+06 0.00118  2.88094E+06 0.00179  3.88152E+06 0.00147  3.57335E+06 0.00123  3.84910E+06 0.00162  2.63446E+06 0.00089  2.30258E+06 0.00176  2.23506E+06 0.00113  1.87535E+06 0.00123  1.78148E+06 0.00160  1.69850E+06 0.00142  1.64259E+06 0.00180  1.55876E+06 0.00128  1.50005E+06 0.00115  1.12477E+06 0.00135  6.91923E+05 0.00173  1.70461E+06 0.00121  1.33235E+06 0.00111  2.53774E+06 0.00130  2.45346E+06 0.00110  1.78926E+06 0.00155  1.16878E+06 0.00119  1.39815E+06 0.00194  1.36991E+06 0.00159  1.18082E+06 0.00147  2.17703E+06 0.00106  4.62339E+05 0.00156  5.73398E+05 0.00154  5.15258E+05 0.00187  2.99928E+05 0.00206  5.16367E+05 0.00150  3.50901E+05 0.00182  3.02737E+05 0.00159  5.87700E+04 0.00119  5.80991E+04 0.00351  5.97172E+04 0.00286  6.11207E+04 0.00205  6.05480E+04 0.00263  5.99054E+04 0.00382  6.18440E+04 0.00320  5.82568E+04 0.00243  1.10358E+05 0.00219  1.78678E+05 0.00108  2.34279E+05 0.00134  7.04146E+05 0.00165  1.05950E+06 0.00133  1.79769E+06 0.00126  1.60423E+06 0.00131  1.33888E+06 0.00172  1.10031E+06 0.00186  1.30093E+06 0.00123  2.40232E+06 0.00162  3.05527E+06 0.00196  5.31333E+06 0.00202  6.94968E+06 0.00167  8.54800E+06 0.00184  4.67395E+06 0.00202  3.05606E+06 0.00181  2.06917E+06 0.00174  1.76826E+06 0.00172  1.66760E+06 0.00188  1.34385E+06 0.00180  8.75194E+05 0.00204  7.88883E+05 0.00185  6.87953E+05 0.00167  5.69852E+05 0.00226  4.33572E+05 0.00169  2.77123E+05 0.00197  9.80088E+04 0.00192 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  2.12617E+17 0.00132  1.98481E+17 0.00178 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.52719E-01 2.0E-05  8.92766E-01 2.7E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  4.90659E-04 0.00038  1.04438E-02 6.5E-05 ];
INF_ABS                   (idx, [1:   4]) = [  4.90659E-04 0.00038  1.04438E-02 6.5E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  8.33512E-08 0.00036  2.23957E-06 6.5E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.52229E-01 2.1E-05  8.82316E-01 2.9E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  2.39021E-01 0.00011  2.77743E-01 4.1E-05 ];
INF_SCATT2                (idx, [1:   4]) = [  9.18806E-02 9.7E-05  1.01232E-01 0.00012 ];
INF_SCATT3                (idx, [1:   4]) = [  4.88562E-03 0.00468  4.13858E-02 0.00050 ];
INF_SCATT4                (idx, [1:   4]) = [ -1.26847E-02 0.00089  1.99566E-02 0.00086 ];
INF_SCATT5                (idx, [1:   4]) = [ -1.86463E-03 0.00377  1.13511E-02 0.00148 ];
INF_SCATT6                (idx, [1:   4]) = [  4.02684E-03 0.00354  7.34902E-03 0.00251 ];
INF_SCATT7                (idx, [1:   4]) = [  7.67013E-04 0.01775  5.23222E-03 0.00198 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.52229E-01 2.1E-05  8.82316E-01 2.9E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.39021E-01 0.00011  2.77743E-01 4.1E-05 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.18806E-02 9.7E-05  1.01232E-01 0.00012 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.88562E-03 0.00468  4.13858E-02 0.00050 ];
INF_SCATTP4               (idx, [1:   4]) = [ -1.26847E-02 0.00089  1.99566E-02 0.00086 ];
INF_SCATTP5               (idx, [1:   4]) = [ -1.86462E-03 0.00377  1.13511E-02 0.00148 ];
INF_SCATTP6               (idx, [1:   4]) = [  4.02683E-03 0.00354  7.34902E-03 0.00251 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.67017E-04 0.01775  5.23222E-03 0.00198 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.29067E-01 0.00019  5.70759E-01 1.5E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.45518E+00 0.00019  5.84017E-01 1.5E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  4.90606E-04 0.00039  1.04438E-02 6.5E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  2.34375E-02 0.00015  1.26399E-02 0.00080 ];

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

INF_CHIT                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.29282E-01 2.6E-05  2.29473E-02 0.00016  2.18972E-03 0.00139  8.80126E-01 3.0E-05 ];
INF_S1                    (idx, [1:   8]) = [  2.31245E-01 0.00012  7.77640E-03 0.00104  9.48814E-04 0.00140  2.76794E-01 4.0E-05 ];
INF_S2                    (idx, [1:   8]) = [  9.27446E-02 5.2E-05 -8.64039E-04 0.00492  3.64145E-04 0.00326  1.00868E-01 0.00012 ];
INF_S3                    (idx, [1:   8]) = [  7.22143E-03 0.00323 -2.33581E-03 0.00081  6.40336E-05 0.01601  4.13218E-02 0.00050 ];
INF_S4                    (idx, [1:   8]) = [ -1.15276E-02 0.00121 -1.15709E-03 0.00389 -4.38074E-05 0.02521  2.00004E-02 0.00083 ];
INF_S5                    (idx, [1:   8]) = [ -1.54283E-03 0.00508 -3.21800E-04 0.00999 -6.30757E-05 0.00597  1.14142E-02 0.00146 ];
INF_S6                    (idx, [1:   8]) = [  4.13482E-03 0.00337 -1.07982E-04 0.02798 -5.43225E-05 0.01427  7.40334E-03 0.00252 ];
INF_S7                    (idx, [1:   8]) = [  8.45258E-04 0.01552 -7.82452E-05 0.00700 -4.03453E-05 0.01708  5.27257E-03 0.00197 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.29282E-01 2.6E-05  2.29473E-02 0.00016  2.18972E-03 0.00139  8.80126E-01 3.0E-05 ];
INF_SP1                   (idx, [1:   8]) = [  2.31245E-01 0.00012  7.77640E-03 0.00104  9.48814E-04 0.00140  2.76794E-01 4.0E-05 ];
INF_SP2                   (idx, [1:   8]) = [  9.27446E-02 5.2E-05 -8.64039E-04 0.00492  3.64145E-04 0.00326  1.00868E-01 0.00012 ];
INF_SP3                   (idx, [1:   8]) = [  7.22143E-03 0.00323 -2.33581E-03 0.00081  6.40336E-05 0.01601  4.13218E-02 0.00050 ];
INF_SP4                   (idx, [1:   8]) = [ -1.15276E-02 0.00121 -1.15709E-03 0.00389 -4.38074E-05 0.02521  2.00004E-02 0.00083 ];
INF_SP5                   (idx, [1:   8]) = [ -1.54282E-03 0.00508 -3.21800E-04 0.00999 -6.30757E-05 0.00597  1.14142E-02 0.00146 ];
INF_SP6                   (idx, [1:   8]) = [  4.13481E-03 0.00337 -1.07982E-04 0.02798 -5.43225E-05 0.01427  7.40334E-03 0.00252 ];
INF_SP7                   (idx, [1:   8]) = [  8.45262E-04 0.01552 -7.82452E-05 0.00700 -4.03453E-05 0.01708  5.27257E-03 0.00197 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  1.13238E-02 0.00140  8.49385E-02 0.00181 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  1.10528E-02 0.00127  8.12213E-02 0.00210 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  1.11486E-02 0.00145  9.06653E-02 0.00200 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  1.17985E-02 0.00152  8.34863E-02 0.00142 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  2.94367E+01 0.00140  3.92446E+00 0.00181 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  3.01584E+01 0.00127  4.10409E+00 0.00210 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  2.98993E+01 0.00144  3.67658E+00 0.00200 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  2.82524E+01 0.00151  3.99270E+00 0.00142 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50813E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19140 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95767E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  2])  = '13' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  2.58604E+05 0.00463  1.21418E+06 0.00217  2.87437E+06 0.00178  3.88597E+06 0.00166  3.58249E+06 0.00215  3.85876E+06 0.00158  2.63869E+06 0.00136  2.30670E+06 0.00127  2.24169E+06 0.00141  1.87933E+06 0.00122  1.78507E+06 0.00142  1.70353E+06 0.00128  1.64693E+06 0.00170  1.56100E+06 0.00158  1.50122E+06 0.00141  1.12814E+06 0.00226  6.92749E+05 0.00151  1.70637E+06 0.00150  1.33292E+06 0.00127  2.54060E+06 0.00157  2.45670E+06 0.00117  1.79168E+06 0.00157  1.17042E+06 0.00164  1.40304E+06 0.00146  1.37514E+06 0.00083  1.18250E+06 0.00127  2.18324E+06 0.00090  4.63569E+05 0.00190  5.74139E+05 0.00114  5.16569E+05 0.00149  2.99884E+05 0.00177  5.18333E+05 0.00190  3.50711E+05 0.00165  3.03135E+05 0.00199  5.88035E+04 0.00184  5.82680E+04 0.00155  5.96910E+04 0.00114  6.15126E+04 0.00124  6.07609E+04 0.00243  6.00956E+04 0.00316  6.18111E+04 0.00300  5.80722E+04 0.00333  1.10776E+05 0.00088  1.79043E+05 0.00316  2.35052E+05 0.00183  7.06570E+05 0.00104  1.06165E+06 0.00145  1.80139E+06 0.00097  1.60655E+06 0.00082  1.34216E+06 0.00086  1.10292E+06 0.00120  1.30180E+06 0.00104  2.40506E+06 0.00082  3.06284E+06 0.00133  5.32617E+06 0.00124  6.96199E+06 0.00116  8.56327E+06 0.00114  4.68382E+06 0.00101  3.06214E+06 0.00122  2.07290E+06 0.00110  1.77199E+06 0.00109  1.67155E+06 0.00110  1.34573E+06 0.00106  8.75885E+05 0.00084  7.90396E+05 0.00102  6.89065E+05 0.00112  5.71507E+05 0.00133  4.34896E+05 0.00137  2.78014E+05 0.00136  9.82028E+04 0.00140 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  2.12993E+17 0.00148  1.98865E+17 0.00114 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.52774E-01 5.1E-05  8.92779E-01 2.0E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  4.90854E-04 0.00037  1.04441E-02 4.7E-05 ];
INF_ABS                   (idx, [1:   4]) = [  4.90854E-04 0.00037  1.04441E-02 4.7E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  8.33960E-08 0.00029  2.23963E-06 4.7E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.52286E-01 4.9E-05  8.82336E-01 1.9E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  2.39048E-01 0.00015  2.77801E-01 0.00011 ];
INF_SCATT2                (idx, [1:   4]) = [  9.18817E-02 0.00021  1.01303E-01 0.00036 ];
INF_SCATT3                (idx, [1:   4]) = [  4.88843E-03 0.00381  4.14226E-02 0.00042 ];
INF_SCATT4                (idx, [1:   4]) = [ -1.27059E-02 0.00155  1.99512E-02 0.00046 ];
INF_SCATT5                (idx, [1:   4]) = [ -1.88755E-03 0.00971  1.13740E-02 0.00144 ];
INF_SCATT6                (idx, [1:   4]) = [  4.00173E-03 0.00235  7.38737E-03 0.00163 ];
INF_SCATT7                (idx, [1:   4]) = [  7.45441E-04 0.01581  5.21297E-03 0.00208 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.52286E-01 4.9E-05  8.82336E-01 1.9E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.39048E-01 0.00015  2.77801E-01 0.00011 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.18817E-02 0.00021  1.01303E-01 0.00036 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.88843E-03 0.00381  4.14226E-02 0.00042 ];
INF_SCATTP4               (idx, [1:   4]) = [ -1.27059E-02 0.00155  1.99512E-02 0.00046 ];
INF_SCATTP5               (idx, [1:   4]) = [ -1.88755E-03 0.00971  1.13740E-02 0.00144 ];
INF_SCATTP6               (idx, [1:   4]) = [  4.00173E-03 0.00235  7.38737E-03 0.00163 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.45441E-04 0.01581  5.21297E-03 0.00208 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.29182E-01 0.00027  5.70723E-01 5.9E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.45445E+00 0.00027  5.84055E-01 5.9E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  4.90833E-04 0.00035  1.04441E-02 4.7E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  2.34572E-02 0.00030  1.26375E-02 0.00057 ];

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

INF_CHIT                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.29317E-01 4.9E-05  2.29690E-02 0.00029  2.19424E-03 0.00094  8.80141E-01 2.1E-05 ];
INF_S1                    (idx, [1:   8]) = [  2.31257E-01 0.00016  7.79102E-03 0.00113  9.48053E-04 0.00311  2.76853E-01 0.00010 ];
INF_S2                    (idx, [1:   8]) = [  9.27392E-02 0.00021 -8.57557E-04 0.00589  3.64364E-04 0.00257  1.00939E-01 0.00036 ];
INF_S3                    (idx, [1:   8]) = [  7.22214E-03 0.00279 -2.33371E-03 0.00107  6.37179E-05 0.00927  4.13589E-02 0.00043 ];
INF_S4                    (idx, [1:   8]) = [ -1.15433E-02 0.00163 -1.16254E-03 0.00106 -4.52559E-05 0.01292  1.99965E-02 0.00047 ];
INF_S5                    (idx, [1:   8]) = [ -1.55897E-03 0.01331 -3.28583E-04 0.01322 -6.46482E-05 0.01777  1.14386E-02 0.00149 ];
INF_S6                    (idx, [1:   8]) = [  4.11007E-03 0.00241 -1.08345E-04 0.01837 -5.42783E-05 0.02291  7.44165E-03 0.00155 ];
INF_S7                    (idx, [1:   8]) = [  8.17256E-04 0.01439 -7.18145E-05 0.03321 -4.04029E-05 0.01517  5.25337E-03 0.00208 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.29317E-01 4.9E-05  2.29690E-02 0.00029  2.19424E-03 0.00094  8.80141E-01 2.1E-05 ];
INF_SP1                   (idx, [1:   8]) = [  2.31257E-01 0.00016  7.79102E-03 0.00113  9.48053E-04 0.00311  2.76853E-01 0.00010 ];
INF_SP2                   (idx, [1:   8]) = [  9.27392E-02 0.00021 -8.57557E-04 0.00589  3.64364E-04 0.00257  1.00939E-01 0.00036 ];
INF_SP3                   (idx, [1:   8]) = [  7.22214E-03 0.00279 -2.33371E-03 0.00107  6.37179E-05 0.00927  4.13589E-02 0.00043 ];
INF_SP4                   (idx, [1:   8]) = [ -1.15433E-02 0.00163 -1.16254E-03 0.00106 -4.52559E-05 0.01292  1.99965E-02 0.00047 ];
INF_SP5                   (idx, [1:   8]) = [ -1.55897E-03 0.01331 -3.28583E-04 0.01322 -6.46482E-05 0.01777  1.14386E-02 0.00149 ];
INF_SP6                   (idx, [1:   8]) = [  4.11008E-03 0.00241 -1.08345E-04 0.01837 -5.42783E-05 0.02291  7.44165E-03 0.00155 ];
INF_SP7                   (idx, [1:   8]) = [  8.17256E-04 0.01439 -7.18145E-05 0.03321 -4.04029E-05 0.01517  5.25337E-03 0.00208 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  1.05938E-02 0.00153  7.24872E-02 0.00111 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  1.03515E-02 0.00143  7.54604E-02 0.00109 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  1.02715E-02 0.00160  6.87594E-02 0.00133 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  1.12078E-02 0.00158  7.35773E-02 0.00096 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  3.14654E+01 0.00153  4.59854E+00 0.00111 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  3.22019E+01 0.00143  4.41735E+00 0.00109 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  3.24527E+01 0.00160  4.84786E+00 0.00133 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  2.97416E+01 0.00158  4.53040E+00 0.00096 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2021 15:03:54' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 30])  = 'hexassembly3D_materialwise.inp' ;
WORKING_DIRECTORY         (idx, [1: 38])  = '/nfs/home/mpt/gen-ffusion/ref-geometry' ;
HOSTNAME                  (idx, [1:  9])  = 'seanode04' ;
CPU_TYPE                  (idx, [1: 47])  = 'AMD Ryzen Threadripper 2990WX 32-Core Processor' ;
CPU_MHZ                   (idx, 1)        = 134251021.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Oct 18 12:54:08 2022' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Oct 18 13:29:13 2022' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1666090448559 ;
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
OMP_THREADS               (idx, 1)        = 50 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  50]) = [  9.65192E-01  9.99677E-01  1.01314E+00  1.01185E+00  1.01590E+00  1.01477E+00  9.98886E-01  9.97976E-01  9.99024E-01  1.00417E+00  1.00341E+00  1.01394E+00  9.99130E-01  9.80655E-01  9.97278E-01  9.86142E-01  9.81509E-01  1.00581E+00  9.98419E-01  9.82590E-01  1.00289E+00  1.00769E+00  9.93051E-01  9.91384E-01  9.93949E-01  1.00609E+00  9.94262E-01  1.00456E+00  9.97045E-01  1.01761E+00  1.00848E+00  1.00167E+00  1.01868E+00  1.00619E+00  9.85597E-01  1.00359E+00  1.01140E+00  9.82246E-01  9.93347E-01  9.95452E-01  1.00175E+00  1.00398E+00  9.84240E-01  1.00182E+00  1.00999E+00  9.86024E-01  1.00790E+00  1.00582E+00  1.01238E+00  1.00144E+00  ];
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

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  3.01306E-01 0.00012  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  6.98694E-01 5.3E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  4.38359E-01 3.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.00693E-01 4.2E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.92888E+00 8.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  5.35330E+01 8.2E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  5.34329E+01 8.3E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  3.55191E+01 0.00010  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.73131E+00 0.00014  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SIMULATED_HISTORIES       (idx, 1)        = 99998525 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  9.99985E+05 0.00018 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.05915E+03 ;
RUNNING_TIME              (idx, 1)        =  3.50813E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.18467E-01  1.18467E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.00000E-03  2.00000E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.49509E+01  3.49509E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.50609E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 30.19140 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.10936E+01 0.00365 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  5.95767E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 128699.88 ;
ALLOC_MEMSIZE             (idx, 1)        = 7733.95;
MEMSIZE                   (idx, 1)        = 7340.32;
XS_MEMSIZE                (idx, 1)        = 505.73;
MAT_MEMSIZE               (idx, 1)        = 27.76;
RES_MEMSIZE               (idx, 1)        = 129.74;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 6677.09;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 393.63;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 32 ;
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

NORM_COEF                 (idx, [1:   4]) = [  7.52110E+10 0.00012  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.65686E-01 0.00036 ];
U235_FISS                 (idx, [1:   4]) = [  2.99633E+16 0.00013  9.71749E-01 2.4E-05 ];
U238_FISS                 (idx, [1:   4]) = [  8.71096E+14 0.00085  2.82507E-02 0.00083 ];
U235_CAPT                 (idx, [1:   4]) = [  8.30040E+15 0.00030  2.25043E-01 0.00029 ];
U238_CAPT                 (idx, [1:   4]) = [  6.34112E+15 0.00033  1.71922E-01 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 99998525 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.91839E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 49015553 4.90402E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40975949 4.09973E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 10007023 1.00117E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 99998525 1.00049E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.78594E-04 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.00000E+06 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.57032E+16 1.5E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.08376E+16 1.2E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  3.68884E+16 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  6.77259E+16 6.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  7.52110E+16 0.00012 ];
TOT_FLUX                  (idx, [1:   2]) = [  7.69020E+18 0.00011 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  7.52991E+15 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  7.52558E+16 8.0E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  4.02068E+18 0.00013 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.95606E+00 9.6E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  5.11008E-01 0.00015 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.48845E-01 9.1E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.03865E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.14756E-01 3.8E-05 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.83741E-01 1.1E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11841E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00644E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45490E+00 1.6E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02399E+02 1.2E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00641E+00 0.00014  9.99403E-01 0.00013  7.03760E-03 0.00165 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00655E+00 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00644E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11840E+00 6.5E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.37192E+01 7.0E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.37185E+01 4.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.20241E-05 0.00095 ];
IMP_EALF                  (idx, [1:   2]) = [  2.20381E-05 0.00066 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.83145E-01 0.00067 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.83181E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.78632E-03 0.00125  2.25808E-04 0.00646  1.19880E-03 0.00321  1.14813E-03 0.00273  2.62180E-03 0.00205  1.12468E-03 0.00330  4.67099E-04 0.00414 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  4.83892E-01 0.00174  1.33445E-02 2.5E-05  3.26735E-02 2.8E-05  1.20922E-01 1.4E-05  3.04281E-01 3.3E-05  8.55735E-01 6.3E-05  2.87446E+00 0.00010 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  7.04113E-03 0.00195  2.32574E-04 0.01073  1.24069E-03 0.00506  1.18767E-03 0.00473  2.72719E-03 0.00330  1.16805E-03 0.00489  4.84948E-04 0.00731 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  4.84363E-01 0.00283  1.33441E-02 3.4E-05  3.26741E-02 4.0E-05  1.20920E-01 2.1E-05  3.04304E-01 5.1E-05  8.55755E-01 9.6E-05  2.87426E+00 0.00015 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  5.11545E-05 0.00044  5.11679E-05 0.00045  4.92541E-05 0.00438 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  5.14824E-05 0.00039  5.14959E-05 0.00040  4.95695E-05 0.00436 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.99325E-03 0.00168  2.32807E-04 0.00993  1.23566E-03 0.00456  1.18235E-03 0.00401  2.69943E-03 0.00274  1.16141E-03 0.00428  4.81600E-04 0.00700 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  4.84217E-01 0.00276  1.33440E-02 3.5E-05  3.26730E-02 4.7E-05  1.20924E-01 2.4E-05  3.04278E-01 4.8E-05  8.55869E-01 0.00010  2.87467E+00 0.00014 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  5.09687E-05 0.00105  5.09806E-05 0.00106  4.93061E-05 0.01242 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.12953E-05 0.00103  5.13074E-05 0.00104  4.96213E-05 0.01240 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  7.06351E-03 0.00563  2.16511E-04 0.03118  1.23920E-03 0.01426  1.16230E-03 0.01394  2.76096E-03 0.01002  1.18033E-03 0.01524  5.04207E-04 0.02147 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.93293E-01 0.00822  1.33432E-02 0.00012  3.26824E-02 0.00014  1.20914E-01 7.6E-05  3.04341E-01 0.00021  8.55735E-01 0.00032  2.87378E+00 0.00051 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  7.06775E-03 0.00566  2.19157E-04 0.02993  1.23543E-03 0.01466  1.16866E-03 0.01349  2.76527E-03 0.00959  1.17728E-03 0.01556  5.01956E-04 0.02094 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.91968E-01 0.00797  1.33431E-02 0.00012  3.26825E-02 0.00014  1.20914E-01 7.3E-05  3.04343E-01 0.00021  8.55681E-01 0.00031  2.87338E+00 0.00050 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.38579E+02 0.00586 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  5.11187E-05 0.00026 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  5.14464E-05 0.00020 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.02637E-03 0.00106 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.37453E+02 0.00112 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.48543E-07 0.00015 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  3.86066E-06 0.00010  3.86032E-06 0.00010  3.91047E-06 0.00146 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.04743E-04 0.00016  1.04753E-04 0.00016  1.03332E-04 0.00166 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.10548E-01 0.00010  5.10554E-01 0.00010  5.09863E-01 0.00228 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07937E+01 0.00289 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  5.34329E+01 8.3E-05  5.26663E+01 0.00013 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  2])  = '14' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  2.56296E+05 0.00056  1.21112E+06 0.00220  2.87571E+06 0.00129  3.88617E+06 0.00139  3.57730E+06 0.00116  3.85592E+06 0.00173  2.63674E+06 0.00183  2.30482E+06 0.00147  2.23864E+06 0.00085  1.87950E+06 0.00091  1.78341E+06 0.00130  1.70058E+06 0.00103  1.64578E+06 0.00123  1.56087E+06 0.00094  1.50116E+06 0.00068  1.12634E+06 0.00100  6.92966E+05 0.00114  1.70877E+06 0.00110  1.33074E+06 0.00077  2.54094E+06 0.00086  2.45574E+06 0.00087  1.78964E+06 0.00074  1.17078E+06 0.00078  1.40080E+06 0.00075  1.37402E+06 0.00064  1.18259E+06 0.00110  2.17990E+06 0.00095  4.62382E+05 0.00055  5.73317E+05 0.00131  5.16139E+05 0.00088  3.00731E+05 0.00108  5.17098E+05 0.00193  3.51394E+05 0.00117  3.02694E+05 0.00053  5.90240E+04 0.00302  5.80704E+04 0.00157  5.95275E+04 0.00247  6.16637E+04 0.00309  6.05617E+04 0.00075  5.99137E+04 0.00077  6.18190E+04 0.00308  5.81983E+04 0.00180  1.10431E+05 0.00137  1.78323E+05 0.00203  2.34621E+05 0.00159  7.07660E+05 0.00146  1.05989E+06 0.00145  1.80127E+06 0.00132  1.60664E+06 0.00146  1.34073E+06 0.00101  1.10253E+06 0.00128  1.30196E+06 0.00099  2.40565E+06 0.00109  3.06090E+06 0.00117  5.32165E+06 0.00082  6.95904E+06 0.00080  8.56067E+06 0.00095  4.68256E+06 0.00083  3.06162E+06 0.00064  2.07194E+06 0.00102  1.77355E+06 0.00065  1.67061E+06 0.00043  1.34631E+06 0.00091  8.76376E+05 0.00099  7.90522E+05 0.00062  6.89261E+05 0.00087  5.70716E+05 0.00037  4.35041E+05 0.00087  2.77694E+05 0.00110  9.81960E+04 0.00119 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  2.12850E+17 0.00102  1.98803E+17 0.00091 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.52827E-01 6.5E-05  8.92800E-01 2.3E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  4.90725E-04 0.00028  1.04447E-02 5.8E-05 ];
INF_ABS                   (idx, [1:   4]) = [  4.90725E-04 0.00028  1.04447E-02 5.8E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  8.33969E-08 0.00025  2.23977E-06 5.8E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.52337E-01 6.3E-05  8.82363E-01 2.9E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  2.39077E-01 0.00015  2.77741E-01 0.00014 ];
INF_SCATT2                (idx, [1:   4]) = [  9.19353E-02 0.00012  1.01233E-01 0.00029 ];
INF_SCATT3                (idx, [1:   4]) = [  4.92946E-03 0.00376  4.14153E-02 0.00063 ];
INF_SCATT4                (idx, [1:   4]) = [ -1.26917E-02 0.00195  1.99637E-02 0.00116 ];
INF_SCATT5                (idx, [1:   4]) = [ -1.88531E-03 0.01342  1.13856E-02 0.00147 ];
INF_SCATT6                (idx, [1:   4]) = [  4.00185E-03 0.00286  7.40275E-03 0.00126 ];
INF_SCATT7                (idx, [1:   4]) = [  7.39466E-04 0.01139  5.21069E-03 0.00213 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.52337E-01 6.3E-05  8.82363E-01 2.9E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.39077E-01 0.00015  2.77741E-01 0.00014 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.19353E-02 0.00012  1.01233E-01 0.00029 ];
INF_SCATTP3               (idx, [1:   4]) = [  4.92946E-03 0.00376  4.14153E-02 0.00063 ];
INF_SCATTP4               (idx, [1:   4]) = [ -1.26917E-02 0.00195  1.99637E-02 0.00116 ];
INF_SCATTP5               (idx, [1:   4]) = [ -1.88531E-03 0.01342  1.13856E-02 0.00147 ];
INF_SCATTP6               (idx, [1:   4]) = [  4.00185E-03 0.00286  7.40275E-03 0.00126 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.39467E-04 0.01139  5.21069E-03 0.00213 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.29191E-01 0.00015  5.70818E-01 6.8E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.45439E+00 0.00015  5.83958E-01 6.8E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  4.90711E-04 0.00028  1.04447E-02 5.8E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  2.34584E-02 6.7E-05  1.26364E-02 0.00082 ];

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

INF_CHIT                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.29369E-01 6.6E-05  2.29689E-02 7.8E-05  2.19911E-03 0.00170  8.80164E-01 3.2E-05 ];
INF_S1                    (idx, [1:   8]) = [  2.31285E-01 0.00017  7.79239E-03 0.00076  9.53305E-04 0.00151  2.76788E-01 0.00014 ];
INF_S2                    (idx, [1:   8]) = [  9.27931E-02 0.00014 -8.57764E-04 0.00600  3.66250E-04 0.00321  1.00867E-01 0.00029 ];
INF_S3                    (idx, [1:   8]) = [  7.26212E-03 0.00279 -2.33265E-03 0.00116  6.15084E-05 0.00804  4.13538E-02 0.00064 ];
INF_S4                    (idx, [1:   8]) = [ -1.15301E-02 0.00200 -1.16155E-03 0.00266 -4.62117E-05 0.02828  2.00099E-02 0.00120 ];
INF_S5                    (idx, [1:   8]) = [ -1.55993E-03 0.01666 -3.25376E-04 0.00936 -6.50776E-05 0.01730  1.14507E-02 0.00143 ];
INF_S6                    (idx, [1:   8]) = [  4.11074E-03 0.00280 -1.08892E-04 0.00922 -5.38829E-05 0.00761  7.45664E-03 0.00129 ];
INF_S7                    (idx, [1:   8]) = [  8.11537E-04 0.01044 -7.20711E-05 0.02985 -4.05397E-05 0.02172  5.25123E-03 0.00199 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.29369E-01 6.6E-05  2.29689E-02 7.8E-05  2.19911E-03 0.00170  8.80164E-01 3.2E-05 ];
INF_SP1                   (idx, [1:   8]) = [  2.31285E-01 0.00017  7.79239E-03 0.00076  9.53305E-04 0.00151  2.76788E-01 0.00014 ];
INF_SP2                   (idx, [1:   8]) = [  9.27931E-02 0.00014 -8.57764E-04 0.00600  3.66250E-04 0.00321  1.00867E-01 0.00029 ];
INF_SP3                   (idx, [1:   8]) = [  7.26211E-03 0.00279 -2.33265E-03 0.00116  6.15084E-05 0.00804  4.13538E-02 0.00064 ];
INF_SP4                   (idx, [1:   8]) = [ -1.15301E-02 0.00200 -1.16155E-03 0.00266 -4.62117E-05 0.02828  2.00099E-02 0.00120 ];
INF_SP5                   (idx, [1:   8]) = [ -1.55993E-03 0.01666 -3.25376E-04 0.00936 -6.50776E-05 0.01730  1.14507E-02 0.00143 ];
INF_SP6                   (idx, [1:   8]) = [  4.11074E-03 0.00280 -1.08892E-04 0.00922 -5.38829E-05 0.00761  7.45664E-03 0.00129 ];
INF_SP7                   (idx, [1:   8]) = [  8.11538E-04 0.01044 -7.20711E-05 0.02985 -4.05397E-05 0.02172  5.25123E-03 0.00199 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  9.93026E-03 0.00097  6.31210E-02 0.00094 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  9.60499E-03 0.00096  6.19379E-02 0.00085 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  9.60669E-03 0.00101  6.19105E-02 0.00106 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  1.06496E-02 0.00096  6.56588E-02 0.00103 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  3.35676E+01 0.00097  5.28088E+00 0.00094 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  3.47043E+01 0.00096  5.38175E+00 0.00085 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  3.46982E+01 0.00100  5.38414E+00 0.00106 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  3.13002E+01 0.00096  5.07677E+00 0.00103 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

