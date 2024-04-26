BEGIN BULK
GRID,1,,0.0
GRID,2,,100.0
GRID,3,,200.0
GRID,4,,300.0
GRID,9,,200.0,200.0
CBEAM,1,30,1,2,9
CBEAM,2,72,2,3,9
CBEAM,3,73,3,4,9
CBEAM,4,31,1,2,9
CBEAM,5,32,1,2,9
$*
$* PROPERTY CARDS
$*
PBEAML        30       2 MSCBML0       L                                +
+       100.0000400.0000 18.0000 13.0000             YES1.000000100.0000+
+       400.0000 18.0000 13.0000
PBEAML        31       2 MSCBML0       L                                +
+       100.0000400.0000 15.0000 15.0000             YES1.000000100.0000+
+       400.0000 15.0000 15.0000
PBEAML        32       2 MSCBML0       L                                +
+       300.0000300.0000 15.0000 15.0000             YES1.000000300.0000+
+       300.0000 15.0000 15.0000
PBEAML        72       2 MSCBML0       T                                +
+       150.0000420.0000 20.0000 11.0000             YES1.000000150.0000+
+       420.0000 20.0000 11.0000
PBEAML        73       2 MSCBML0       I                                +
+       150.0000 50.0000 60.0000 10.0000 15.0000 20.0000             YES+
+       1.000000150.0000 50.0000 60.0000 10.0000 15.0000 20.0000
$*
$* MATERIAL CARDS
$*
MAT1           22.0600+87.9231+70.3000007.8500-6
ENDDATA