%chk=75-89-8.chk
%mem=4GB
%nprocshared=4
# opt freq b3lyp/aug-cc-pvtz EmpiricalDispersion=GD3 geom=connectivity

Title Card Required

0 1
 C                 -0.39291963   -0.06933062   -0.04403881
 C                  0.45214225   -0.31707329   -1.25863505
 H                 -1.26196936   -0.77079538   -0.01729347
 H                  0.20319935   -0.21680190    0.88917972
 O                 -0.89868825    1.26824109   -0.04403881
 H                 -1.43453970    1.40822104   -0.82817454
 F                  1.50204999    0.53096038   -1.29076967
 F                  0.92992967   -1.57969715   -1.25906965
 F                 -0.26860897   -0.13866058   -2.38610415
 O                 -2.59346053    1.31238894   -2.29200680
 H                 -2.13401928    0.99764820   -3.07396077
 H                 -3.43785637    1.69117884   -2.54718101

 1 2 1.0 3 1.0 4 1.0 5 1.0
 2 7 1.0 8 1.0 9 1.0
 3
 4
 5 6 1.0
 6
 7
 8
 9
 10 11 1.0 12 1.0
 11
 12

