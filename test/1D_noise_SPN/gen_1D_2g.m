
f = 1.0; % Hz
mkdir input


%% Geometry data
DZ = 5.0;  % cm
DX=NaN;
DY=NaN;

L = 300.0; % cm
z_prima = -32.5;
n_nodes = floor(L / DZ);

save input/GEOM_data DX DY DZ

%% Cross Section data

% Fuel 1
D1_F     = 1.0/(3.0 * 4.67810e-01);
D2_F     = 1.0/(3.0 * 1.09620e-01);
ABS1_F   = 1.17659e-02;
ABS2_F   = 1.07186e-01;
NUFIS1_F = 5.62285e-03;
NUFIS2_F = 1.45865e-01;
REM_F    = 1.60795e-02;

D1     = D1_F     * ones(1,1, n_nodes);
D2     = D2_F     * ones(1,1, n_nodes);
ABS1   = ABS1_F   * ones(1,1, n_nodes);
ABS2   = ABS2_F   * ones(1,1, n_nodes);
NUFIS1 = NUFIS1_F * ones(1,1, n_nodes);
NUFIS2 = NUFIS2_F * ones(1,1, n_nodes);
REM    = REM_F    * ones(1,1, n_nodes);

save input/XS_data D1 D2 ABS1 ABS2 REM NUFIS1 NUFIS2

%% Precursors data

v1 = 1.25e7;
v2 = 2.5e5;
Beff=0.0065;
l = 0.0784;

save input/DYN_data Beff l v1 v2 f


%% Set pertubation data

dABS1 = zeros(1,1,n_nodes);
dABS2= zeros(1,1,n_nodes);
dREM = zeros(1,1,n_nodes);
dNUFIS1 = zeros(1,1,n_nodes);
dNUFIS2= zeros(1,1,n_nodes);

cell_prima = floor((L/2+z_prima)/DZ) + 1;
dABS1(cell_prima) = ABS1_F/DZ/1000;
dABS2(cell_prima) = ABS2_F/DZ/1000;

save input/dS_data dNUFIS1 dNUFIS2 dABS1 dABS2 dREM z_prima
