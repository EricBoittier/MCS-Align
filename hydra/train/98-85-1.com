%chk=98-85-1.chk
%mem=5GB
%nprocshared=4
# opt freq b3lyp/aug-cc-pvtz EmpiricalDispersion=GD3 geom=connectivity

Title Card Required

0 1
 C                  2.15872300    0.38924300   -1.04583700
 C                  0.78072500    0.43330500   -1.23939200
 C                 -0.09261600    0.37362400   -0.15711400
 C                  0.43636000    0.25910000    1.12909600
 C                  1.80919500    0.21463300    1.32585600
 C                  2.67570900    0.28022200    0.23801800
 H                  2.82317600    0.43389600   -1.89552100
 H                  0.37849000    0.50901000   -2.24012800
 H                 -0.24211200    0.18618700    1.96540000
 H                  2.20633800    0.12654100    2.32607800
 H                  3.74333200    0.24215500    0.39160800
 C                 -1.59203800    0.41587100   -0.36163500
 O                 -2.25683500   -0.54841600    0.43511400
 H                 -1.87597000   -1.41044600    0.21593900
 C                 -2.16699700    1.77115400    0.01003900
 H                 -1.96989600    1.97920400    1.05992900
 H                 -1.71731000    2.55758100   -0.59237600
 H                 -3.24287100    1.76920200   -0.14867700
 O                 -0.52911900   -2.72689000   -0.31549700
 H                 -0.23051400   -3.17423800    0.47851500
 H                  0.15081600   -2.06837600   -0.48969300
 H                 -1.79021200    0.22881900   -1.42220100

 1 2 1.5 6 1.5 7 1.0
 2 3 1.5 8 1.0
 3 4 1.5 12 1.0
 4 5 1.5 9 1.0
 5 6 1.5 10 1.0
 6 11 1.0
 7
 8
 9
 10
 11
 12 13 1.0 15 1.0 22 1.0
 13 14 1.0
 14
 15 16 1.0 17 1.0 18 1.0
 16
 17
 18
 19 20 1.0 21 1.0
 20
 21
 22

