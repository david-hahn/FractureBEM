rem we use a low-frequency layered material to bias seeding towards the middle and propagation towards the base of the cube
..\bin\FractureBEM.exe cube 3[0,0,-5e-3] 2[0,0,5e-3] ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 2e6 -C 2  ^
  -M lay(0,0,1,1,1.570796326,5e5/0,0,1,1,0,76e5) ^
  -s 30 -r 0.1 -R 0.005 -o _results/cube_bias -x 1 -v 2

pause
rem last tested 21-Jul-2015