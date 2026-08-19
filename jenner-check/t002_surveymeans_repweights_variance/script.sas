/* Adapted from TonyAnSGF2020_sample.sas — the sample script's final analysis
   step (`proc surveymeans ... varmethod=bootstrap; repweights RepWt:; run;`)
   run against the repo's own Mortality survey extract, supplying replicate
   weight columns directly via REPWEIGHTS= (the macro's own user-supplied-
   replicate-weights mode, surveycalibrate.sas lines 302-309) rather than
   through the macro's internal OUTWEIGHTS= generation step. */

data Mortality;
   input ID VarStrata VarPSU SWeight Age VitalStatus PovArInd Gender
         RepWt_1 RepWt_2 RepWt_3 RepWt_4;
   datalines;
      1  03  1  13312    66  1   1   1   12800   13900   13000   13600
      2  03  1   7941    71  3   1   2    7600    8300    7800    8100
      3  03  1  16048     .  4   1   1   15400   16700   15700   16400
      4  03  3   9298    58  3   1   1    8900    9700    9100    9500
      5  03  2  15336    56  3   1   2   14700   15900   15000   15700
      6  03  1  14744    63  1   1   1   14100   15300   14400   15100
      7  03  2  83729    70  1   2   2   80100   87400   82200   85300
      8  03  3 106492    57  1   2   1  102000  111100  104600  108500
      9  03  3  78083    81  3   2   2   74800   81600   76700   79600
     10  03  3  55957    79  3   2   1   53600   58500   55000   57000
     11  03  3  83729    68  1   2   2   80100   87400   82200   85300
     12  03  1  78083    78  3   2   2   74800   81600   76700   79600
     13  03  2  13824    78  1   1   2   13200   14400   13500   14100
     14  03  3  13824    70  3   1   2   13200   14400   13500   14100
     15  03  3  44649    50  1   1   2   42800   46700   44000   45500
     16  03  1   9298     .  6   1   1    8900    9700    9100    9500
     17  03  1  13824    77  1   1   2   13200   14400   13500   14100
     18  03  3   4767    82  3   1   1    4600    5000    4700    4900
     19  03  3  15336    56  3   1   2   14700   15900   15000   15700
     20  03  3  16048    68  3   1   1   15400   16700   15700   16400
;

proc print data=Mortality (obs=10); run;

proc surveymeans data = Mortality mean;
   weight SWeight;
   class VitalStatus;
   var VitalStatus Age;
   repweights RepWt_1-RepWt_4;
run;
