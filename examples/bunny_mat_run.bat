
..\bin\FractureBEM.exe ^
  bunny 1(fixed) 2[0.0005,-0.0023,0] ^
  --remesh 1000 -O 3 --max-cracks 6 ^
  -s 20 -r 0.010 -R 0.0005 ^
  -E 3100e6 -n 0.327 -d 1200 ^
  -S 76e6 -K 5e5 --compressive 2 ^
  -o _results\bunny_mat --no-si ^
  -M vdb(mat_grains.vdb,4e5,1e-4,1)

pause

