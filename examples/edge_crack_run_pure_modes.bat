@echo off
rem edge crack test case
rem want to show pure mode fracture surfaces (I, II, III)

rem pure mode I case
..\bin\FractureBEM.exe edge_crack 1[0,0,-5e-4] 2[0,0,5e-4] 5(crack) -c 6 ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 1e6 --compressive 2 -x 0 ^
  -s 10 -r 0.15 -R 0.003 -o _results\edg_crk_modeI

rem pure mode I case compressive (since the material is tougher under compression, this should not crack completely)
..\bin\FractureBEM.exe edge_crack 1[0,0,5e-4] 2[0,0,-5e-4] 5(crack) -c 6 ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 1e6 --compressive 2 -x 0 ^
  -s 10 -r 0.15 -R 0.003 -o _results\edg_crk_modeIc

rem pure mode II case
..\bin\FractureBEM.exe edge_crack 1(fixed) 2(0,5e6,0) 5(crack) -c 6 ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 1e6 --compressive 2 -x 0 ^
  -s 10 -r 0.15 -R 0.003 -o _results\edg_crk_modeII

rem pure mode III case
..\bin\FractureBEM.exe edge_crack 3(8e6,0,0) 4(-8e6,0,0) 5(crack) -c 6 ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 1e6 --compressive 2 -x 0 ^
  -s 10 -r 0.15 -R 0.003 -o _results\edg_crk_modeIII


pause
rem last tested 21-Jul-2015