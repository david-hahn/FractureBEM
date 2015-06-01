@echo off
rem edge crack test case
rem want to show pure mode fracture surfaces (I, II, III) (note that III needs mode I contribution to avoid compression)

rem pure mode I case
..\bin\FractureBEM.exe edge_crack 1[0,0,-5e-4] 2[0,0,5e-4] 5(crack) -c 6 ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 1e6 --compressive 2 -x 0 ^
  -s 10 -r 0.15 -R 0.003 -o _results\edg_crk_modeI

rem pure mode I case compressive
rem ..\bin\FractureBEM.exe edge_crack 1[0,0,5e-4] 2[0,0,-5e-4] 5(crack) -c 6 ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 1e6 --compressive 2 -x 0 ^
  -s 10 -r 0.15 -R 0.003 -o _results\edg_crk_modeIc

rem pure mode II case
..\bin\FractureBEM.exe edge_crack 1(fixed) 2(0,5e6,0) 5(crack) -c 6 ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 1e6 --compressive 2 -x 0 ^
  -s 10 -r 0.15 -R 0.003 -o _results\edg_crk_modeII

rem (almost) pure mode III case -- mode I contribution required to avoid compression on the crack
..\bin\FractureBEM.exe edge_crack 3(5e6,0,1e4) 4(-5e6,0,-1e4) 5(crack) 9(fixed) -c 6 ^
  -E 3100e6 -n 0.327 -d 1200 -S 76e6 -K 1e6 --compressive 2 -x 0 ^
  -s 10 -r 0.15 -R 0.003 -o _results\edg_crk_modeIII


pause
