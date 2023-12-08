#$ -cwd              # Ejecutar en directorio actual.
#$ -l h_vmem=100g     # Memoria requerida.
#$ -l h_rt=48:00:00  # Tiempo requerido (hh:mm:ss).

echo '3D_Seaborg'
./femffusion.exe -f 3D_Seaborg/periferal_fast/periferal_seaborg_SP1_p10.prm
./femffusion.exe -f 3D_Seaborg/periferal_fast/periferal_seaborg_SP1_p5.prm
./femffusion.exe -f 3D_Seaborg/seaborg_SP1_fe1_p20.prm
./femffusion.exe -f 3D_Seaborg/seaborg_SP1_fe1_p10.prm
./femffusion.exe -f 3D_Seaborg/seaborg_SP1_fe1_p40.prm