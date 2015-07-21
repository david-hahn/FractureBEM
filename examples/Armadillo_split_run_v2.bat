..\bin\FractureBEM.exe Armadillo.ply 1(fixed) 2(-5.5e5,0,0) 3(5.5e5,0,0) ^
  --remesh 1000 -x 0 ^
  -M vdb(mat_rolled_grains_1.vdb,1e5,0.1,1) ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 1e6 -C 2 ^
  -s 30 -r 10 -R 0.2 -o _results\arm_split ^
  -X (0,70,0,1,0,0,0,-1,0) -v 2 --no-si


pause
rem last tested 21-Jul-2015: 10 steps, segments with threshold 0.1
