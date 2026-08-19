/* Adapted from TonyAnSGF2020_sample.sas — the repo's own Mortality survey
   extract (from the public SGF 2020 paper) run through PROC SURVEYMEANS,
   the design-based summary step the sample script performs before handing
   off to %SurveyCalibrate. Exercises stratified/clustered weighted survey
   estimation exactly as authored, unmodified. */

data Mortality;
   input ID VarStrata VarPSU SWeight Age VitalStatus PovArInd Gender;
   datalines;
      1  03  1  13312    66  1   1   1
      2  03  1   7941    71  3   1   2
      3  03  1  16048     .  4   1   1
      4  03  3   9298    58  3   1   1
      5  03  2  15336    56  3   1   2
      6  03  1  14744    63  1   1   1
      7  03  2  83729    70  1   2   2
      8  03  3 106492    57  1   2   1
      9  03  3  78083    81  3   2   2
     10  03  3  55957    79  3   2   1
     11  03  3  83729    68  1   2   2
     12  03  1  78083    78  3   2   2
     13  03  2  13824    78  1   1   2
     14  03  3  13824    70  3   1   2
     15  03  3  44649    50  1   1   2
     16  03  1   9298     .  6   1   1
     17  03  1  13824    77  1   1   2
     18  03  3   4767    82  3   1   1
     19  03  3  15336    56  3   1   2
     20  03  3  16048    68  3   1   1
     21  03  1   9298    74  1   1   1
     22  03  2  14744     .  6   1   1
     23  03  2   4767    77  3   1   1
     24  03  2  16048    65  3   1   1
     25  03  1 106492    61  1   2   1
     26  03  3 170748     .  1   2   2
     27  03  2   9298     .  1   1   1
     28  03  1  78083    89  1   2   2
     29  03  1 170748    58  1   2   2
     30  03  2  20029    64  1   1   2
;

proc print data=Mortality (obs=10); run;

proc surveymeans data = Mortality sum;
   class PovArInd Gender;
   weight SWeight;
   var PovArInd Gender;
run;

proc surveymeans data = Mortality mean;
   weight SWeight;
   class VitalStatus;
   var VitalStatus Age;
   cluster VarPSU;
   strata VarStrata;
run;
