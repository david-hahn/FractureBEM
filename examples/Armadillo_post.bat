..\bin\VisMesh.exe _results\arm_lay -o _view\arm_lay -s 6 -S 15 ^
  -r 10 -q 0.001 -Q 1e-11 --vis-obj --vis-close

..\bin\VisMesh.exe _results\arm_gr -o _view\arm_gr -s 13 ^
  -r 10 -q 0.001 -Q 1e-11 -P 1e-12 --vis-obj --vis-close --segment-last 0.1 --segment-files

..\bin\VisMesh.exe _results\arm_split -o _view\arm_split -s 10 ^
  -r 10 -q 0.001 -Q 1e-11 -P 1e-12 --vis-obj --vis-close --segment-last 0.1 --segment-files

pause
