
..\bin\FractureBEM.exe notched_bar 1(0,2e7,0) 2(fixed) 3(fixed) ^
  -E 3100e6 -n 0.327 -d 1200 ^
  -S 76e6 -K 1e6 --compressive 2 ^
  -M vdb(mat_grains.vdb,1e5,3e-5,0) -x 5 ^
  -s 20 -r 4e-3 -R 2e-4 -o _results\bar_grains

..\bin\FractureBEM.exe notched_bar 1(0,2e7,0) 2(fixed) 3(fixed) ^
  -E 3100e6 -n 0.327 -d 1200 ^
  -S 76e6 -K 1e6 --compressive 2 ^
  -M vdb(mat_grains.vdb,1e5,3e-5,0) -x 5 -v 2^
  -s 20 -r 4e-3 -R 2e-4 -o _results\bar_gr_v2

pause