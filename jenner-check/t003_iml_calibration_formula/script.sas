/* Excerpted from surveycalibrate.sas's %NONE calibration branch (lines
   337-389): the core generalized-least-squares calibration weight formula
   -- CalW = w + w#(X*beta), beta = ginv(X'*wx)*(T - X'*w) -- run standalone
   against the repo's own Mortality survey extract (from the public SGF 2020
   paper) with the Poverty/NonPoverty/Male/Female control totals the sample
   script uses. The macro's internal bookkeeping (CALL SYMPUTX writeback to
   communicate status back to the wrapping %macro) is omitted here since
   this bundle calls the IML calibration math directly rather than through
   the full %SurveyCalibrate macro. */

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
;

data Mortality; set Mortality;
   _intercept_=1;
   Poverty=0; NonPoverty=0; Male=0; Female=0;
   if (Gender=1)   then Male      =1;
   if (Gender=2)   then Female    =1;
   if (PovArInd=1) then Poverty   =1;
   if (PovArInd=2) then NonPoverty=1;
run;

proc iml;
  use Mortality;
  read all var {SWeight} into w;
  read all var {_intercept_ Poverty NonPoverty Male Female} into x;
  close Mortality;
  T={95000 60000 35000 40000 55000}`;
  wx=w#x;
  beta=ginv(X`*wx)*(T-X`*w);
  CalW=w+w#(X*beta);
  print beta[label="Calibration coefficients (beta)"];
  print CalW[label="Calibrated weights"];
  create CalWeights from CalW[colname="Cal_SWeight"];
  append from CalW;
  quit;

proc print data=CalWeights; run;
