%chk=80-73-9.chk
%mem=4GB
%nprocshared=4
# opt freq b3lyp/aug-cc-pvtz EmpiricalDispersion=GD3 geom=connectivity

Title Card Required

0 1
 C                 -0.18679323    1.33055292   -0.05095489
 C                  1.25210887    1.36931805   -0.55371339
 C                  0.86239683    3.14735682    0.82582375
 H                 -0.88732452    0.55262453   -0.27227667
 H                  0.17452020    0.77398900    0.78844324
 H                  1.26914869    1.87040373   -1.49897673
 H                  1.70185212    0.40569968   -0.67233318
 N                 -0.45234609    2.46058512    0.84301995
 N                  1.92544405    2.18545330    0.46503480
 O                  1.04592183    4.36503729    1.08494423
 C                  3.16265792    2.78278903   -0.05783407
 H                  2.93644268    3.37384482   -0.92060850
 H                  3.60449338    3.40320399    0.69367555
 H                  3.84759474    2.00611366   -0.32716088
 C                 -1.57247456    3.29976202    0.39358214
 H                 -1.37498801    3.65876570   -0.59485876
 H                 -2.47249922    2.72110802    0.38968879
 H                 -1.68526806    4.13024178    1.05877459
 O                 -0.74111281    5.93110022    1.54772500
 H                  0.17584264    5.66110050    1.63656681
 H                 -0.79129562    6.88922075    1.51475937

 1 2 1.0 4 1.0 5 1.0 8 1.0
 2 6 1.0 7 1.0 9 1.0
 3 8 1.0 9 1.0 10 2.0
 4
 5
 6
 7
 8 15 1.0
 9 11 1.0
 10
 11 12 1.0 13 1.0 14 1.0
 12
 13
 14
 15 16 1.0 17 1.0 18 1.0
 16
 17
 18
 19 20 1.0 21 1.0
 20
 21

