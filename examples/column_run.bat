..\bin\FractureBEM.exe Column.stl --remesh 1000 -O 5 1(fixed) 2[0.1,-0.05,0] ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 1e4 -C 2 -x 17 ^
  -s 30 -r 0.1 -R 0.0086 -o _results\column ^
  -M vdb(mat_grains.vdb,0.2e4,0.001,0)


pause