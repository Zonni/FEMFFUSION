#$ -cwd              # Ejecutar en directorio actual.

echo '3D_Seaborg central p5'
./femffusion.exe -f 3D_Seaborg/seaborg_SP1_p5.prm

echo '3D_Seaborg central p10'
./femffusion.exe -f 3D_Seaborg/seaborg_SP1_p10.prm

echo '3D_Seaborg periferal fast p20'
./femffusion.exe -f 3D_Seaborg/periferal_fast/periferal_seaborg_SP1_p20.prm

echo '3D_Seaborg periferal slow p20'
./femffusion.exe -f 3D_Seaborg/periferal_slow/seaborg_SP1_p20.prm

echo '3D_Seaborg periferal fast p10'
./femffusion.exe -f 3D_Seaborg/periferal_fast/periferal_seaborg_SP1_p10.prm

echo '3D_Seaborg periferal slow p10'
./femffusion.exe -f 3D_Seaborg/periferal_slow/seaborg_SP1_p10.prm
