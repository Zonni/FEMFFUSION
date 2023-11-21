#!/bin/bash
set -e  # Stop at first error

echo 'REACTOR 1D_hom_slab of 2cm'
echo '   1D_hom_zerocurrent'
./femffusion.exe -f 1D/1D_hom_vacuum.prm
echo '   1D_hom_zerocurrent'
./femffusion.exe -f 1D/1D_hom_zerocurrent.prm
echo '   1D_hom_zeroflux'
./femffusion.exe -f 1D/1D_hom_zeroflux.prm
echo ''
echo ''
echo 'REACTOR 2D_biblis'
./femffusion.exe -f 2D_BIBLIS/biblis_SP1.prm
./femffusion.exe -f 2D_BIBLIS/biblis_SP3.prm
./femffusion.exe -f 2D_BIBLIS/biblis_SP5.prm
echo ''
echo ''
echo 'REACTOR 2D_C5G7'
./femffusion.exe -f 2D_C5G7/c5g7_ref1_FE2.prm -solver_type bifpam
./femffusion.exe -f 2D_C5G7/c5g7_SP3_ref1_FE2.prm
echo ''
echo ''
echo 'REACTOR 2D_VVER440'
./femffusion.exe -f 2D_VVER440/2D_VVER440_FE3.prm
./femffusion.exe -f 2D_VVER440/2D_VVER440_SP3_FE3.prm
echo ''
echo ''
echo 'REACTOR 3D_IAEA'
./femffusion.exe -f 3D_IAEA/IAEA_FE3.prm
./femffusion.exe -f 3D_IAEA/IAEA_SP3_FE3.prm
echo ''
echo ''
echo 'REACTOR 3D_LANGENBUCH'
./femffusion.exe -f 3D_Langenbuch/3D_Langenbuch_fe3.prm
echo ''
echo ''
echo 'REACTOR 3D_VVER440'
./femffusion.exe -f 3D_VVER440/3D_VVER440_FE3.prm
