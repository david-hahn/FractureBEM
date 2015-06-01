..\bin\FractureBEM.exe Armadillo.ply 1(fixed) 2(-5e5,0,0) 3(5e5,0,0) ^
  --remesh 1000 --max-cracks 6 ^
  -M vdb(mat_grains.vdb,9e5,0.05,1) ^
  -E 3100e6 -n 0.327 -d 1200 ^
  -S 76e6 -K 1e6 --compressive 2 ^
  -s 12 -r 3 -R 0.2 -o _results\arm_gr 

..\bin\FractureBEM.exe Armadillo.ply 1(fixed) 2(-5e5,0,0) 3(5e5,0,0) ^
  --remesh 1000 ^
  -M lay(0.2,1,0,1,0,1.6e6) ^
  -E 3100e6 -n 0.327 -d 1200 ^
  -S 76e6 -K 1e6 --compressive 2 ^
  -s 10 -r 5 -R 0.3 -o _results\arm_lay --out-sub

pause