%chk=62-53-3.chk
%mem=2GB
%nprocshared=2
# opt freq b3lyp/aug-cc-pvtz EmpiricalDispersion=GD3 geom=connectivity

Title Card Required

0 1
 C                  2.21566900   -0.73488100   -0.20322600
 C                  1.14973900   -1.40874300    0.38236300
 C                 -0.04702000   -0.75542100    0.63686000
 C                 -0.18878200    0.59517900    0.31513300
 C                  0.87803700    1.27053500   -0.27475700
 C                  2.06975100    0.60668600   -0.53215400
 H                  3.14396300   -1.24749500   -0.40027000
 H                  1.24686700   -2.45262900    0.64009500
 H                 -0.88506900   -1.28769200    1.06166600
 H                  0.77244200    2.31571100   -0.52968500
 H                  2.88758600    1.14449100   -0.98748800
 N                 -1.42390300    1.23660600    0.52030200
 H                 -1.89919100    0.92794500    1.35376600
 H                 -1.36290100    2.24150900    0.50349200
 O                 -3.14438200   -0.80749400   -0.67867900
 H                 -2.58134900   -1.29358100   -1.28183600
 H                 -2.66432800    0.01531800   -0.51773800

 1 2 1.5 6 1.5 7 1.0
 2 3 1.5 8 1.0
 3 4 1.5 9 1.0
 4 5 1.5 12 1.0
 5 6 1.5 10 1.0
 6 11 1.0
 7
 8
 9
 10
 11
 12 13 1.0 14 1.0
 13
 14
 15 16 1.0 17 1.0
 16
 17

