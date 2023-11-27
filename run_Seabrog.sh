#$ -cwd              # Ejecutar en directorio actual.
#$ -l h_vmem=100g     # Memoria requerida.
#$ -l h_rt=48:00:00  # Tiempo requerido (hh:mm:ss).

echo '3D_Seaborg'
./femffusion.exe -f 3D_Seaborg/seaborg_SP1_fe2.prm -type_prec bad-broyden -init_prec gs-cgilu
