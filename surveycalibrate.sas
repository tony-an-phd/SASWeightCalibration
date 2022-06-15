
/*--------------------------------------------------------------------

Author: Tony An, SAS Institute Inc. 
        tony.an@sas.com

In survey sampling, calibration is commonly used for adjusting weights
to ensure that estimates for covariates in the sample match known
auxiliary information, such as marginal totals from census
data. Calibration can also be used to adjust for unit
nonresponse. This macro for calibration was developed using SAS/STAT
15.1 and SAS/IML software. The macro enables you to input the design
information, the controls for the auxiliary variables, and your
preferred calibration method, including either linear or exponential
method. Because unbounded calibration method can result in extreme
calibration weights, this macro also supports bounded versions of both
linear and exponential calibration methods. The macro creates
calibration replication weights according to the sample design and the
specified calibration method.

----------------------------------------------------------------------
 DISCLAIMER:

       THIS INFORMATION IS PROVIDED BY SAS INSTITUTE INC. AS A SERVICE
 TO ITS USERS.  IT IS PROVIDED "AS IS".  THERE ARE NO WARRANTIES,
 EXPRESSED OR IMPLIED, AS TO MERCHANTABILITY OR FITNESS FOR A
 PARTICULAR PURPOSE REGARDING THE ACCURACY OF THE MATERIALS OR CODE
 CONTAINED HEREIN.

---------------------------------------------------------------------*/

%macro SurveyCalibrate(
   DATA=,      /* Input data set name                                        */
   OUT=,       /* Output data set name                                       */
   /*  Calibration parameters                                                */
   METHOD=,    /* LINEAR | EXPONENTIAL | TRUNLINEAR | LOGIT                  */
   WEIGHT=,    /* Original weight variable                                   */
   CALWT=,     /* Calibration weight variable, default CalWt                 */
   CONTROLVAR=,/* Auxiliary control variables for calibration                */
   CTRLTOTAL=, /* Marginal totals for CONTROLVAR                             */
   EPS=,       /* Convergence criterion for stopping iteration, default=0.01 */
   MAXITER=,   /* Maximum number of iteration, default=25                    */
   LOWER=,     /* Lower bound, must be in (0,1)                              */
   UPPER=,     /* Upper bound, must be bigger than 1 or .                    */
   NOINT=,     /* Do not keep sum of sampling weights unchanged              */
   /* Replication parameters                                                 */
   NOREPWT=,   /* Request no replicate weights                               */
   VARMETHOD=, /* BRR | JK | BOOTSTRAP, default is JK                        */
   REPS=,      /* Number of replicates for bootstrap or brr                  */
   CLUSTER=,   /* Cluster variables                                          */
   STRATA=,    /* Strata variables                                           */
   SEED=,      /* Random seed                                                */
   FAY=,       /* Fay coefficient for BRR varmethod                          */
   RATE=,      /* FPC for bootstrap replicate weights                        */
   OUTJKCOEFS=,/* OUTJKCOEFS data set                                        */                                               
   REPWEIGHTS= /* Replicate weight variables                                 */        
   );


   /* ********** initilization *********** */
   %let switchToExp=0;
   %let negtiveCalW=0;
   %let lowerBSearchNeeded=0;
   %let upperBSearchNeeded=0;  
   %let maxIterReached=0;
   %let success=1;
   %let ceateIntercept=1;
   %let creatRepWt=1;
   %let calRepWt=1;
   %let badrep=0;
   %let lowestbound=0.001;

   ods graphics off;
   
   %if &data eq %then %do;
      %put ERROR: No SAS data set is specified.;
      %let success=0;
      %goto EXIT;
   %end;
  
   %if &out eq %then %do;
      %put ERROR: You must provide an out= data set name.;
      %let success=0;
      %goto EXIT;
   %end;

   %if &weight eq %then %do;
      %put ERROR: No weight variable is specified.;
      %let success=0;
      %goto EXIT;
   %end;
  
   %if &controlvar eq %then %do;
      %put ERROR: You must provide CONTROLVAR= variables.;
      %let success=0;
      %goto EXIT;
   %end;
  
   %if &ctrltotal eq %then %do;
      %put ERROR: You must provide CTRLTOTAL= totals.;
      %let success=0;
      %goto EXIT;
   %end;
  
   %let nctrlvar=%sysfunc(countw(&controlvar, %str( )));
   %let nctr=%sysfunc(countw(&ctrltotal, %str( )));
 
   %if &nctrlvar ne &nctr %then %do;
      %put ERROR: The number of CONTROLVAR= variables &nctrlvar and the number of CTRLTOTAL= totals &nctr do not match.;
      %let success=0;
      %goto EXIT;
   %end; 
  
   %if &calwt eq %then %let calwt=Cal_&weight;
  
   %if &eps eq %then %let eps=0.01;
   %if ((&eps<=0) | &eps>=1) %then %do;
     %put ERROR: The eps=&eps is outside the range of (0,1).;
     %let success=0;
     %goto EXIT;
   %end;
  
   %if &method eq %then %do;
      %let method=NONE;
   %end;
   %else %do;
      %let method=%upcase(&method); 
      %if ((&method^=LINEAR) & (&method^=EXPONENTIAL) &
       (&method^=TRUNLINEAR) & (&method^=LOGIT) 
       & (&method^=NONE)) %then %do;
         %put ERROR: method=&method is invalid.;
         %let success=0;
         %goto EXIT; 
         %end;
      %end;
      
   %if &norepwt ne %then %do;
      %if ((%upcase(&norepwt)^=TRUE) & (%upcase(&norepwt)^=FALSE)) %then %do;
         %put ERROR: NOREPWT=&norepwt is invalid macro parameter, it only takes keywords TRUE or FALSE.;
         %let creatRepWt=0;
         %let success=0;   
         %goto EXIT;
         %end;
      %if (%upcase(&norepwt)=TRUE) %then %do;
         %let calRepWt=0;
         %let creatRepWt=0;
         %end;
      %end;

   %if &maxiter ne %then %do;
      %if (&maxiter<2) %then %do;
         %put ERROR: MAXITER=&maxiter is too small.;
         %let success=0;   
         %goto EXIT;
         %end;
      %end;
   %if &maxiter eq %then %let maxiter=50;      

   %let lowerPrivided=0;
   %if &lower eq %then %let lower=&lowestbound;
      %else %do;
         %let lowerPrivided=1;
         %if ((&lower>=1) | (&lower<=0)) %then %do;
            %put ERROR: The lower=&lower is outside the range of (0,1).;
            %let success=0;
            %goto EXIT;
         %end;
      %end;
   
   %if (&lowerPrivided=0 & ((&method=TRUNLINEAR)|(&method=LOGIT)))
     %then %let lowerBSearchNeeded=1;
   
   %let upperPrivided=0;
   %if &upper eq %then %let upper=.;
      %else %do;
         %let upperPrivided=1;
         %if (&upper<1) %then %do;
            %if (&upper ^=. ) %then %do;
               %put ERROR: The upper=&upper is smaller than 1.;
               %let success=0; 
               %goto EXIT;
            %end;
         %end;
      %end;
   %if (&upperPrivided=0 & ((&method=TRUNLINEAR)|(&method=LOGIT)))
      %then %do;
      %let upperBSearchNeeded=1;
      /* get initial value for the upperbound from data */
      ods listing close; 
      proc surveymeans data=&data min max; 
         var &weight; 
         ods output statistics=_weightminmax_;
      run;
      data _weightminmax_; 
         set _weightminmax_; 
         _maxminratio_=max/min;
         call symputx('upper', _maxminratio_);
      run;
      %end; 

   %if ((&method=NONE)|(&method=LINEAR)|(&method=EXPONENTIAL))
     %then %do;
       %let lower=1e-4;
       %let upper=.;
     %end;

   %if (%upcase(&noint)=TRUE) %then %let ceateIntercept=0;

   %if ((&calRepWt=1)&(&creatRepWt=1)) %then %do;
       %if ((&varmethod ne ) & (&repweights eq ))
          %then %let varmethod=%upcase(&varmethod);
       %if (&varmethod eq ) %then %let varmethod=JK;
       %if ((&varmethod^=BRR) & (&varmethod^=JK) &
           (&varmethod^=JACKKNIFE) &
           (&varmethod^=BOOTSTRAP)) %then %do;
          %put ERROR: varmethod=&varmethod is invalid.;
          %let creatRepWt=0;
          %let success=0;
          %goto EXIT;
          %end;
       %end;
   
   /* remove missing values, creat intercept */   
   data _obs_&data; set &data;
      _obs_=_n_;
   data _temp_&data; set _obs_&data;
      %do i=1 %to &nctrlvar;
         if (missing(%scan(&controlvar, &i))=1) then delete;
      %end;
      %if (&ceateIntercept=1) %then %do;
         _intercept_=1;
         %end; 
    run;
   /* get the total weights under clean data */
   %if (&ceateIntercept=1) %then %do;
      ods listing close;
      proc surveymeans data=_temp_&data sumwgt; 
         var _intercept_; 
         weight &weight; 
         ods output statistics=_temp_sumwgt_;
      run;
      data _temp_sumwgt_; set _temp_sumwgt_; 
         call symputx('totalSumWt', sumwgt);
      run;
      /* add intercept to the control */
      %let controlvar=%str(_intercept_ &controlvar);
      %let ctrltotal=%str(&totalSumWt &ctrltotal);
      /* %put contrvars=&controlvar;
         %put ctrltotal=&ctrltotal; */
   %end; 

   /* for replication method */
   %if (&calRepWt=1) %then %do;
      %let repwtPrefix=RepWt_;
      %if &repweights ne %then %do;
         %let creatRepWt=0;
         %let repwtPrefix=&repweights;
         %end;
      %end;
        
    /* Create replicate weights  */
    ods listing close;
    %if (&creatRepWt=1) %then %do;
       %if ((&varmethod=JK)|(&varmethod=JACKKNIFE)) %then %do;
          %if (&outjkcoefs ne ) %then %let outjkcoefs=%str(outjkcoefs=&outjkcoefs);
          proc surveymeans data=_temp_&data varmethod=jk(&outjkcoefs outweights=_temp_outrepwts);
          %end;
       %if (&varmethod=BRR) %then %do;
          %if (&fay ne ) %then %let fay=%str(fay=&fay);
          proc surveymeans data=_temp_&data varmethod=brr(&fay outweights=_temp_outrepwts);
          %end;
       %if (&varmethod=BOOTSTRAP) %then %do;
          %if (&rate ne ) %then %let rate=%str(r=&rate);
          %if (&seed ne ) %then %let seed=%str(seed=&seed);
          %if (&reps ne ) %then %let reps=%str(reps=&reps);
          proc surveymeans data=_temp_&data &rate
          varmethod=bootstrap(&reps &seed outweights=_temp_outrepwts);
          %end;          
       %if &strata ne %then %do;
          strata &strata;
          %end;
       %if &cluster ne %then %do;
          cluster &cluster;
          %end;
          weight &weight;   
          var _obs_;
          ods output VarianceEstimation=_temp_varestsummary; 
       run;
       %if ((&syserr ne 0) & (&syserr ne 4)) %then %do;
          %let success=0;
          %goto EXIT;
          %end;
       data _temp_varestsummary; set _temp_varestsummary; 
            if (label1 ^= 'Number of Replicates') then delete; 
                call symputx('nReps', cValue1);
       run;
     ods listing close; 
     
     data _temp_&data; merge _temp_&data _temp_outrepwts;
     run;
     %end;
     /* user provided replicate weights */
     %if ((&creatRepWt=0) & (&calRepWt=1)) %then %do;
        proc surveymeans data=_temp_&data varmethod=&varmethod;
        var _obs_;
        repweights &repweights:;
        weight &weight;
        ods output VarianceEstimation=_temp_varestsummary; 
        run;
        %if ((&syserr ne 0) & (&syserr ne 4)) %then %do;
           %let success=0;
           %goto EXIT;
           %end;
        data _temp_varestsummary; set _temp_varestsummary; 
            if (label1 ^= 'Number of Replicates') then delete; 
                call symputx('nReps', cValue1);
        run;
        %end;

OPTIONS errors=0 MERGENOBY=NOWARN ;

   %if ((&method=TRUNLINEAR)|(&method=LOGIT))
      %then %do;
         %if ((&upperBSearchNeeded=1) & (&lowerBSearchNeeded=1))
            %then %goto LU_&method;
         %if ((&upperBSearchNeeded=1) & (&lowerBSearchNeeded=0))
            %then %goto U_&method;
         %if ((&upperBSearchNeeded=0) & (&lowerBSearchNeeded=1))
            %then %goto L_&method;
         %if ((&upperBSearchNeeded=0) & (&lowerBSearchNeeded=0))
            %then %goto &method;
      %end;
      %else %goto &method;

/* method= is not specified */
%NONE:
proc iml;
  switchToExp=0;
  negtiveCalW=0;
  use _temp_&data;;
  read all var {&weight} into w;
  read all var {&controlVar} into x;
  close _temp_&data; 
  T={&ctrltotal}`; 
  wx=w#x; 
  beta=ginv(X`*wx)*(T-X`*w);
  CalW=w+w#(X*beta);
  minCalWValue=min(CalW);
  *print minCalWValue;
  if (minCalWValue<0)then
    negtiveCalW=1;
    start expodist(v) global(w);
       weps = 1e-4;
       verybig= weps**(log(weps)-1);
       f=0;
       do i = 1 to nrow(w);
       if (w > weps) then
          vdw = v[i]/w[i];
       else 
          vdw = v[i]/weps;
       if (vdw > weps) then
          fi = 1 + vdw**(log(vdw) - 1);
       else 
          fi = verybig;  
       f=f+fi;
       end; 
       /* y = (J(1, nrow(w), 1))+ (v/w`)##((log(v/w`)-J(1, nrow(w), 1))); */
       return(f);
    finish expodist;
    /* if have negative weights, use exponential */
    if (negtiveCalW=1) then 
      do;
        switchToExp=1;
        j=J(nrow(T), 1, 0) ;
        constrain=x`||j|| T;
        w0 = w;
        optn = {0 0};
        blc  = constrain;
        l=&lower; u=.;
        lbound=l*w`|| {. .};  
        ubound=u*w`|| {. .};
        bound=lbound//ubound; 
        blc=bound//constrain;
        call nlpqn(rc, v, "expodist", w0, optn, blc);
        Calw=v`;  
        minCalWValue=min(CalW);
      end;
  create _tmpnewwt_ from CalW[colname="&calwt"]; 
  append from CalW;
  call symputx('negtiveCalW', negtiveCalW);
  call symputx('switchToExp', switchToExp);
  quit;
  %goto CALREPWTS_LINEAR;
/* end NONE method */

/* LINEAR method */
%LINEAR:
proc iml;
   negtiveCalW=0;
   use _temp_&data;;
   read all var {&weight} into w;
   read all var {&controlVar} into x;
   close _temp_&data; 
   T={&ctrltotal}`; 
   wx=w#x; 
   beta=ginv(X`*wx)*(T-X`*w);
   CalW=w+w#(X*beta);
   minCalWValue=min(CalW);
   if (minCalWValue<0)then
     negtiveCalW=1;
   create _tmpnewwt_ from CalW[colname="&calwt"]; 
   append from CalW;
   call symputx('negtiveCalW', negtiveCalW);        
   quit;
   %goto CALREPWTS_LINEAR;
/* end of LINEAR method */

/* TRUNLINEAR method */
%TRUNLINEAR:
   proc iml;
   use _temp_&data;; read all var {&weight} into w;
   *print w;
   read all var {&controlVar} into x;
   *print x;
   close _temp_&data; 
   T={&ctrltotal}`;   
   *print T;
   j=J(nrow(T), 1, 0) ;
   constrain=x`||j|| T;

   w0 = w;
   optn = {0 0};
   blc  = constrain;

   start dist(v) global(w);
      y = (v/w`-J(1, nrow(w), 1))`;
      f=sum(w#(y##2));
      return(f);
   finish dist;
 
   start fitwithbound(ll, uu) global(w, optn, constrain);
      lbound=ll*w`|| {. .};  
      ubound=uu*w`|| {. .};
      bound=lbound//ubound;    
      blc=bound//constrain;
      call nlpqn(rc, v, "dist", w, optn, blc);
      return (rc);
   finish fitwithbound;

   *multiplier=(sqrt(5)-1)/2;
   multiplier=0.5;
   ll=&lower; uu=&upper;
   *print ll uu ;
   /*  initial */      
   rc=fitwithbound(ll, uu); 

   /* if specified both bounds, make it or break it */
   if (rc<0) then do;
      success=0;
      call symputx('success', success);
   end;
   else do;
      lbound=ll*w`|| {. .};  
      ubound=uu*w`|| {. .};
      bound=lbound//ubound;  
      blc=bound//constrain;
      call nlpqn(rc, v, "dist", w, optn, blc);
      Calw=v`;  
      create _tmpnewwt_ from CalW[colname="&calwt"]; 
      append from CalW;    
   end;

  %goto CALREPWTS_TRUNLINEAR;
/* end of truncatelinear */
   
/*U_TRUNLINEAR method needs to search upper bound */
%U_TRUNLINEAR:
   proc iml;
   use _temp_&data;; read all var {&weight} into w;
   *print w;
   read all var {&controlVar} into x;
   *print x;
   close _temp_&data;
   T={&ctrltotal}`;   
   *print T;
   j=J(nrow(T), 1, 0) ;
   constrain=x`||j|| T;

   maxIterReached=0;
   w0 = w;
   optn = {0 0};
   blc  = constrain;

   start dist(v) global(w);
      y = (v/w`-J(1, nrow(w), 1))`;
      f=sum(w#(y##2));
      return(f);
   finish dist;
 
   start fitwithbound(ll, uu) global(w, optn, constrain);
      lbound=ll*w`|| {. .};  
      ubound=uu*w`|| {. .};
      bound=lbound//ubound;    
      blc=bound//constrain;
      call nlpqn(rc, v, "dist", w, optn, blc);
      return (rc);
   finish fitwithbound;

   *multiplier=(sqrt(5)-1)/2;
   multiplier=0.5;
   ll=&lower; uu=&upper;
   *print "initial ll and uu=" ll uu ;
   /*  initial */      
   rc=fitwithbound(ll, uu); 

   /* only lower bound is specified, search the upper one */
   success=1;
   **** search for upper bound ****;
   **** default bound is l=0 u=max/min ****;
   initialL=&lower; initialU=&upper; finalrc=0;
   lastuu=uu;
   lastrc=rc;
   absolutebound=1;

   if (rc<=0) then 
     do;
       uu=(1+multiplier)*lastuu; 
       absolutebound=absolutebound<>lastuu;
     end;

   if (rc>0) then
     uu=multiplier*lastuu<>absolutebound; 

   continue=1;
   iter=0;
   do while (continue=1);  
      rc=fitwithbound(ll, uu); 
      if (rc<=0) then absolutebound=absolutebound<>uu;
      temp=uu;
      max=uu<>lastuu;
      min=uu><lastuu;
      mim=min<>absolutebound;  
      diff=max-min;   
      if (rc<0 & lastrc<0) then
        uu=(1+multiplier)*max;
      /* *(max-absolutebound)*multiplier;  */
      if (rc>0 & lastrc>0) then
        uu=min-multiplier*(min-absolutebound); 
      if (rc>0 & lastrc<0) then
        uu=min+(1-multiplier)*diff;
      if (rc<0 & lastrc>0) then 
        uu=min+multiplier*diff; 
      /* always go above absolutebound */
      uu=uu<>absolutebound;   

      lastuu=temp;
      lastrc=rc;

      if (diff<&eps & rc>0) then 
        do;    
          uu=lastuu;
          continue=0;
        end;
      if (diff=0 & rc<0) then do;
        ll=lastll;
        continue=0;
        success=0;
      end;
      iter=iter+1;
      *print iter rc lastrc ll uu lastuu min max absolutebound diff;
      if (iter>&maxiter & continue=1) then do;
        continue=0;
        maxIterReached=1;
        success=0;
      end;
   end; /* end of do while (continue=1); for looking for upper bound */
   if (success=1) then do;
     *print "The lower bound=" ll "Upper bound=" uu;    
     lbound=ll*w`|| {. .};  
     ubound=uu*w`|| {. .};
     bound=lbound//ubound;    
     blc=bound//constrain;
     call nlpqn(rc, v, "dist", w0, optn, blc);
     calw=v`;        
     create _tmpnewwt_ from CalW[colname="&calwt"]; 
     append from CalW;  
   end;
   call symputx('lower', ll);
   call symputx('upper', uu);
   call symputx('printiter', iter);
   call symputx('maxIterReached', maxIterReached);
   call symputx('success', success);
   
  %goto CALREPWTS_TRUNLINEAR;
  /* end of U_TRUNLINEAR */

/* L_TRUNLINEAR method needs to search lower bound */
%L_TRUNLINEAR:
   proc iml;
   use _temp_&data;; read all var {&weight} into w;
   *print w;
   read all var {&controlVar} into x;
   *print x;
   close _temp_&data; 
   T={&ctrltotal}`;   
   *print T;
   j=J(nrow(T), 1, 0) ;
   constrain=x`||j|| T;

   maxIterReached=0;
   w0 = w;
   optn = {0 0};
   blc  = constrain;

   start dist(v) global(w);
      y = (v/w`-J(1, nrow(w), 1))`;
      f=sum(w#(y##2));
      return(f);
   finish dist;
 
   start fitwithbound(ll, uu) global(w, optn, constrain);
      lbound=ll*w`|| {. .};  
      ubound=uu*w`|| {. .};
      bound=lbound//ubound;    
      blc=bound//constrain;
      call nlpqn(rc, v, "dist", w, optn, blc);
      return (rc);
   finish fitwithbound;

   *multiplier=(sqrt(5)-1)/2;
   multiplier=0.5;
   ll=&lower; uu=&upper;
   *print ll uu ;
   /*  initial */      
   rc=fitwithbound(ll, uu); 

   /* only upper bound is specified, search for a lower one */
   success=1;
   **** search for lower bound ****;
   **** default bound is &lowestbound u=max/min ****;
   initialL=&lower; initialU=&upper; finalrc=0;
   lastll=ll;
   lastrc=rc;
   absolutebound=1;
   max=ll<>lastll;
   min=ll><lastll;
   /* since the initial l&lowestbound if LOWER= is not specified, 
      if we cannot get a solution, then no need to go any 
      further, just quit */ 
   if (rc<=0) then success=0;

   if (rc>0) then ll=lastll+multiplier*(absolutebound-lastll);  
  
   continue=1;
   iter=0; 
   do while (continue=1 & success=1);  
     *print "before the next iteration......";
     *print iter rc lastrc min max uu absolutebound;    
     rc=fitwithbound(ll, uu); 
     if (rc<=0) then absolutebound=absolutebound><ll;
     temp=ll;
     max=ll<>lastll;
     min=ll><lastll;
     *print lastll ll min max;
     diff=max-min; 
     if (rc<0 & lastrc<0) then 
       ll=&lowestbound+multiplier*(min-&lowestbound); 
     if (rc>0 & lastrc>0) then 
       ll=max+multiplier*(absolutebound-max); 
     if (rc>0 & lastrc<0) then 
       ll=min+multiplier*diff; 
     if (rc<0 & lastrc>0) then 
       ll=min+(1-multiplier)*diff;
     /* always go under absolutebound */
     ll=ll><absolutebound;

     *print "new rc ll min max diff" rc ll min max diff;
     lastll=temp;
     lastrc=rc;
     *print "going to the next loop" rc lastrc lastll ll;
     
     if (diff<&eps & rc>0) then do;
       ll=lastll;
       continue=0;
     end;
     if (diff=0 & rc<0) then do;
       ll=lastll;
       continue=0;
       success=0;
     end;
     iter=iter+1;
     if (iter>&maxiter & continue=1) then do;
       continue=0;
       maxIterReached=1;
       success=0;
     end;
   end; /* end while (continue=1); searching for lower bound */
   if (success=1) then do;
     *print "The lower bound=" ll "Upper bound=" uu;    
     lbound=ll*w`|| {. .};  
     ubound=uu*w`|| {. .};
     bound=lbound//ubound;    
     blc=bound//constrain;
     call nlpqn(rc, v, "dist", w0, optn, blc);
     calw=v`;        
     create _tmpnewwt_ from CalW[colname="&calwt"]; 
     append from CalW;  
   end;
   call symputx('lower', ll);
   call symputx('upper', uu);
   call symputx('printiter', iter);
   call symputx('maxIterReached', maxIterReached);
   call symputx('success', success);
   
   %goto CALREPWTS_TRUNLINEAR;
   /* end of L_TRUNLINEAR */


/* need to search both upper and lower bounds */
%LU_TRUNLINEAR:
   proc iml;
   use _temp_&data;; read all var {&weight} into w;
   *print w;
   read all var {&controlVar} into x;
   *print x;
   close _temp_&data; 
   T={&ctrltotal}`;   
   *print T;
   j=J(nrow(T), 1, 0) ;
   constrain=x`||j|| T;

   w0 = w;
   optn = {0 0};
   blc  = constrain;
   maxIterReached=0;

   start dist(v) global(w);
      y = (v/w`-J(1, nrow(w), 1))`;
      f=sum(w#(y##2));
      return(f);
   finish dist;
 
   start fitwithbound(ll, uu) global(w, optn, constrain);
      lbound=ll*w`|| {. .};  
      ubound=uu*w`|| {. .};
      bound=lbound//ubound;    
      blc=bound//constrain;
      call nlpqn(rc, v, "dist", w, optn, blc);
      return (rc);
   finish fitwithbound;

   *multiplier=(sqrt(5)-1)/2;
   multiplier=0.5;
   ll=&lower; uu=&upper;
   *print ll uu ;
   /*  initial */      
   rc=fitwithbound(ll, uu); 

   success=1;
   /* search the upper bound first */
       **** search for upper bound ****;
       **** default bound is l=0 u=max/min ****;
       initialL=&lower; initialU=&upper; finalrc=0;
       lastuu=uu;
       lastrc=rc;
       absolutebound=1;

       if (rc<=0) then 
         do;
           uu=(1+multiplier)*lastuu; 
           absolutebound=absolutebound<>lastuu;
         end;

       if (rc>0) then  
         do;
           uu=multiplier*lastuu; 
         end;
       continue=1;
       iter=0;
       do while (continue=1);  
         rc=fitwithbound(ll, uu); 
         if (rc<=0) then do;
           absolutebound=absolutebound<>uu;
         end;
         temp=uu;
         max=uu<>lastuu;
         min=uu><lastuu;
         diff=max-min;   
         if (rc<0 & lastrc<0) then 
           uu=(1+multiplier)*max;
         /* max+(max-absolutebound)*multiplier; */ 
         if (rc>0 & lastrc>0) then 
           uu=min-multiplier*(min-absolutebound); 
         if (rc>0 & lastrc<0) then 
           uu=min+(1-multiplier)*diff; 
         if (rc<0 & lastrc>0) then 
           uu=min+multiplier*diff; 
         /* always go above absolutebound */
         uu=uu<>absolutebound;
      
         lastuu=temp;
         lastrc=rc;
         *print iter rc lastrc ll uu diff;

         if (abs(diff)<&eps & rc>0) then 
           do;    
             uu=lastuu;
             continue=0;
          end;
         iter=iter+1;
         *print continue iter;
         if (iter>&maxiter & continue=1) then do;
           continue=0;
           maxIterReached=1;
           success=0;
         end;
       end; /* end of do while (continue=1); looking for upper bound */

       /* now the upper bound has been set for default lower=0 */

    /* now search for a lower bound */
    printiterl=iter;
    if (success=1) then do;
      **** search for lower bound ****; 
      lastll=ll;
      lastrc=rc;
      absolutebound=1;
      ll=lastll+multiplier*(absolutebound-lastll);  /* lastll should be 0.01 and rc>0 */
      continue=1;
      do while (continue=1);      
        rc=fitwithbound(ll, uu); 
        if (rc<=0) then absolutebound=absolutebound><ll;
        temp=ll;
        max=ll<>lastll;
        min=ll><lastll;
        diff=max-min; 
        if (rc<0 & lastrc<0) then 
          ll=&lowestbound+multiplier*(min-&lowestbound); 
        if (rc>0 & lastrc>0) then 
          ll=max+multiplier*(absolutebound-max); 
        if (rc>0 & lastrc<0) then 
          ll=min+multiplier*diff; 
        if (rc<0 & lastrc>0) then 
          ll=min+(1-multiplier)*diff; 
        /* always go under absolutebound */
        ll=ll><absolutebound;

        lastll=temp;
        lastrc=rc;
        
        *print iter rc lastrc ll uu diff;
        
        if (diff<&eps & rc>0) then do;
          ll=lastll;
          continue=0;
        end;
        if (diff=0 & rc<0) then do;
          ll=lastll;
          continue=0;
          success=0;
        end;
     
        iter=iter+1;
        if (iter>&maxiter+printiterl & continue=1) then do;
          /* since lower=0 is okay, this almost never happen */
          continue=0;
          maxIterReached=1;
          success=0;
        end;    
      end; /* end while (continue=1); searching for lower bound */
    end; /* end of search lower end */
    
  if (success=1) then do;
     *print "The lower bound=" ll "Upper bound=" uu;    
     lbound=ll*w`|| {. .};  
     ubound=uu*w`|| {. .};
     bound=lbound//ubound;    
     blc=bound//constrain;
     call nlpqn(rc, v, "dist", w0, optn, blc);
     calw=v`;        
     create _tmpnewwt_ from CalW[colname="&calwt"]; 
     append from CalW;  
   end;
   call symputx('lower', ll);
   call symputx('upper', uu);
   call symputx('printiter', iter);
   call symputx('maxIterReached', maxIterReached);
   call symputx('success', success);

   %goto CALREPWTS_TRUNLINEAR;
/* end of LU_TRUNLINEAR */

/* EXPONENTIAL method */
%EXPONENTIAL:
proc iml;
   use _temp_&data;; 
   read all var {&weight} into w;
   read all var {&controlVar} into x;
   close _temp_&data; 
   T={&ctrltotal}`; 
   j=J(nrow(T), 1, 0) ;
   constrain=x`||j|| T;
   w0 = w;
   optn = {0 0};
   blc  = constrain;
   
   start dist(v) global(w);
      weps = 1e-4;
      verybig= weps**(log(weps)-1);
      f=0;
      do i = 1 to nrow(w);
      if (w > weps) then
         vdw = v[i]/w[i];
      else 
         vdw = v[i]/weps;
      if (vdw > weps) then
         fi = 1 + vdw**(log(vdw) - 1);
      else 
         fi = verybig;  
      f=f+fi;
      end; 
/* y = (J(1, nrow(w), 1))+ (v/w`)##((log(v/w`)-J(1, nrow(w), 1))); 
      *print y;
      f=sum(y); */
      return(f);
   finish dist;
 
   l=&lower; u=&upper;
   lbound=l*w`|| {. .};  
   ubound=u*w`|| {. .};
   bound=lbound//ubound; 
   blc=bound//constrain;
   call nlpqn(rc, v, "dist", w0, optn, blc);
   Calw=v`;  
   create _tmpnewwt_ from CalW[colname="&calwt"]; 
   append from CalW;
%goto CALREPWTS_LOGIT;  
/* end of EXPONENTIAL method */

/* LOGIT method */
%LOGIT:
proc iml;
   use _temp_&data;; 
   read all var {&weight} into w;
   read all var {&controlVar} into x;
   close _temp_&data; 
   T={&ctrltotal}`; 
   j=J(nrow(T), 1, 0) ;
   constrain=x`||j|| T;
   optn = {0 0};
   blc  = constrain;
   
   start dist(v) global(w);
      weps = 1e-4;
      verybig= weps**(log(weps)-1);
      f=0;
      do i = 1 to nrow(w);
      if (w > weps) then
         vdw = v[i]/w[i];
      else 
         vdw = v[i]/weps;
      if (vdw > weps) then
         fi = 1 + vdw**(log(vdw) - 1);
      else 
         fi = verybig;  
      f=f+fi;
      end;
      /*y = (J(1, nrow(w), 1))+ (v/w`)##((log(v/w`)-J(1, nrow(w), 1)));
      *print y;
      f=sum(y); */
      return(f);
   finish dist;

   w0=w;
   success=1; 
   ll=&lower; uu=&upper;
   lbound=ll*w`|| {. .};  
   ubound=uu*w`|| {. .};  
   bound=lbound//ubound;    
   blc=bound//constrain;
   call nlpqn(rc, v, "dist", w0, optn, blc);

   *print "for LOGIT method" ll uu rc;    
   if (rc<0) then success=0; 
   if (success=1) then do;
     Calw=v`;  
     create _tmpnewwt_ from CalW[colname="&calwt"]; 
     append from CalW;
   end;
   call symputx('lower', ll);
   call symputx('upper', uu);
   call symputx('success', success);
%goto CALREPWTS_LOGIT;  
/* end of LOGIT method */


/* L_LOGIT method need to search lower bound */
%L_LOGIT:
proc iml;
   use _temp_&data;; 
   read all var {&weight} into w;
   read all var {&controlVar} into x;
   close _temp_&data; 
   T={&ctrltotal}`; 
   j=J(nrow(T), 1, 0) ;
   constrain=x`||j|| T;
   optn = {0 0};
   blc  = constrain;
   w0=w;
   
   start dist(v) global(w);
      weps = 1e-4;
      verybig= weps**(log(weps)-1);
      f=0;
      do i = 1 to nrow(w);
      if (w > weps) then
         vdw = v[i]/w[i];
      else 
         vdw = v[i]/weps;
      if (vdw > weps) then
         fi = 1 + vdw**(log(vdw) - 1);
      else 
         fi = verybig;  
      f=f+fi;
      end;
      /* y = (J(1, nrow(w), 1))+ (v/w`)##((log(v/w`)-J(1, nrow(w), 1)));
      *print y;
      f=sum(y); */
      return(f);
   finish dist;

   start fitwithbound(ll, uu) global(w, optn, constrain);
      lbound=ll*w`|| {. .};  
      ubound=uu*w`|| {. .};
      bound=lbound//ubound;    
      blc=bound//constrain;
      call nlpqn(rc, v, "dist", w, optn, blc);
      return (rc);
   finish fitwithbound;

   success=1;
   maxIterReached=0; 
   multiplier=0.5;
   ll=&lower; uu=&upper;
   /*  initial */      
   rc=fitwithbound(ll, uu); 
   initialL=&lower; initialU=&upper; finalrc=0;
   lastll=ll;
   lastrc=rc;
   absolutebound=1;
   /* since the initial ll=&lowestbound if LOWER= is not specified, 
      if we cannot get a solution, then no need to go any 
      further, just quit */ 
   if (rc<=0) then success=0; 
   if (rc>0) then ll=lastll+multiplier*(absolutebound-lastll);  
   continue=1;
   iter=0; 
   max=ll<>lastll;
   min=ll><lastll;
    
   *print "going to search lower bound";
   *print ll uu iter rc;
   do while (continue=1 & success=1);      
     *print iter success ll uu optn blc w;
     rc=fitwithbound(ll, uu); 
     if (rc<=0) then absolutebound=absolutebound><ll;
     temp=ll;
     max=ll<>lastll;
     min=ll><lastll;
     diff=max-min; 
     if (rc<0 & lastrc<0) then 
       ll=&lowestbound+multiplier*(min-&lowestbound); 
     if (rc>0 & lastrc>0) then 
       ll=max+multiplier*(absolutebound-max); 
     if (rc>0 & lastrc<0) then 
       ll=min+multiplier*diff; 
     if (rc<0 & lastrc>0) then 
       ll=min+(1-multiplier)*diff;
     /* always go under absolutebound */
     ll=ll><absolutebound;
     lastll=temp;
     lastrc=rc;
     *print iter rc lastrc ll uu lastll min max absolutebound diff;
     if (diff<&eps & rc>0) then do;
       ll=lastll;
       continue=0;
     end;
     if (diff=0 & rc<0) then do;
       ll=lastll;
       continue=0;
       success=0;
     end;
     iter=iter+1;
     if (iter>&maxiter & continue=1) then do;
       continue=0;
       maxIterReached=1;
       success=0;
     end;
   end; /* end while (continue=1); searching for lower bound */
   if (success=1) then do;
     *print "The lower bound=" ll "Upper bound=" uu;    
     lbound=ll*w`|| {. .};  
     ubound=uu*w`|| {. .};
     bound=lbound//ubound;    
     blc=bound//constrain;
     call nlpqn(rc, v, "dist", w0, optn, blc);
     Calw=v`;  
     create _tmpnewwt_ from CalW[colname="&calwt"]; 
     append from CalW;
   end;
   call symputx('lower', ll);
   call symputx('upper', uu);
   call symputx('printiter', iter);
   call symputx('maxIterReached', maxIterReached);
   call symputx('success', success);
%goto CALREPWTS_LOGIT;  
/* end of L_LOGIT method */

/* U_LOGIT method need to search lower bound */
%U_LOGIT:
proc iml;
   use _temp_&data;; 
   read all var {&weight} into w;
   read all var {&controlVar} into x;
   close _temp_&data; 
   T={&ctrltotal}`; 
   j=J(nrow(T), 1, 0) ;
   constrain=x`||j|| T;
   optn = {0 0};
   blc  = constrain;
   w0=w;
   
   start dist(v) global(w);
      weps = 1e-4;
      verybig= weps**(log(weps)-1);
      f=0;
      do i = 1 to nrow(w);
      if (w > weps) then
         vdw = v[i]/w[i];
      else 
         vdw = v[i]/weps;
      if (vdw > weps) then
         fi = 1 + vdw**(log(vdw) - 1);
      else 
         fi = verybig;  
      f=f+fi;
      end;
      /* y = (J(1, nrow(w), 1))+ (v/w`)##((log(v/w`)-J(1, nrow(w), 1)));
      *print y;
      f=sum(y); */
      return(f);
   finish dist;

   start fitwithbound(ll, uu) global(w, optn, constrain);
      lbound=ll*w`|| {. .};  
      ubound=uu*w`|| {. .};
      bound=lbound//ubound;    
      blc=bound//constrain;
      *print optn blc;
      call nlpqn(rc, v, "dist", w, optn, blc);
      return (rc);
   finish fitwithbound;

   success=1;
   maxIterReached=0; 
   multiplier=0.5;
   ll=&lower; uu=&upper;
   /*  initial */      
   rc=fitwithbound(ll, uu); 
   initialL=&lower; initialU=&upper; finalrc=0;
   lastuu=uu;
   lastrc=rc;
   absolutebound=1;

   if (rc<=0) then do;
     uu=(1+multiplier)*lastuu; 
     absolutebound=absolutebound<>lastuu;
   end;

   if (rc>0) then
     uu=multiplier*lastuu<>absolutebound; 

   continue=1;
   iter=0;
   max=uu<>lastuu;
   min=uu><lastuu;
   diff=max-min;   
   do while (continue=1);      
      *print iter success ll uu optn blc w;
      *print iter rc lastrc ll uu lastuu min max absolutebound diff;
      rc=fitwithbound(ll, uu);
      if (rc<=0) then absolutebound=absolutebound<>uu;
      temp=uu;
      max=uu<>lastuu;
      min=uu><lastuu;
      diff=max-min;     

      if (rc<0 & lastrc<0) then
        uu=(1+multiplier)*max;
      /* *(max-absolutebound)*multiplier;  */
      if (rc>0 & lastrc>0) then
        uu=min-multiplier*(min-absolutebound); 
      if (rc>0 & lastrc<0) then
        uu=min+(1-multiplier)*diff;
      if (rc<0 & lastrc>0) then 
        uu=min+multiplier*diff; 
      /* always go above absolutebound */
      uu=uu<>absolutebound;  

      lastuu=temp;
      lastrc=rc;

      if (diff<&eps & rc>0) then 
        do;    
          uu=lastuu;
          continue=0;
        end;
      if (diff=0 & rc<0) then do;
       ll=lastll;
       continue=0;
       success=0;
      end;
      iter=iter+1;
      *print iter rc lastrc ll uu lastuu min max absolutebound diff;
      if (iter>&maxiter & continue=1) then do;
        continue=0;
        maxIterReached=1;
        success=0;
      end;
   end; /* end of do while (continue=1); for looking for upper bound */
   if (success=1) then do;
     *print "The lower bound=" ll "Upper bound=" uu;    
     lbound=ll*w`|| {. .};  
     ubound=uu*w`|| {. .};
     bound=lbound//ubound;    
     blc=bound//constrain;
     call nlpqn(rc, v, "dist", w0, optn, blc);
     Calw=v`;  
     create _tmpnewwt_ from CalW[colname="&calwt"]; 
     append from CalW;
   end;
   call symputx('lower', ll);
   call symputx('upper', uu);
   call symputx('printiter', iter);
   call symputx('maxIterReached', maxIterReached);
   call symputx('success', success);
%goto CALREPWTS_LOGIT;  
/* end of U_LOGIT method */

/* need to search both upper and lower bounds */
%LU_LOGIT:
proc iml;
   use _temp_&data;; 
   read all var {&weight} into w;
   read all var {&controlVar} into x;
   close _temp_&data; 
   T={&ctrltotal}`; 
   j=J(nrow(T), 1, 0) ;
   constrain=x`||j|| T;
   optn = {0 0};
   blc  = constrain;
   w0 = w;
   start dist(v) global(w);
      weps = 1e-4;
      verybig= weps**(log(weps)-1);
      f=0;
      do i = 1 to nrow(w);
      if (w > weps) then
         vdw = v[i]/w[i];
      else 
         vdw = v[i]/weps;
      if (vdw > weps) then
         fi = 1 + vdw**(log(vdw) - 1);
      else 
         fi = verybig;  
      f=f+fi;
      end;    
      /*y = (J(1, nrow(w), 1))+ (v/w`)##((log(v/w`)-J(1, nrow(w), 1)));
      *print y;
      f=sum(y); */
      return(f);
   finish dist;

   start fitwithbound(ll, uu) global(w, optn, constrain);
      lbound=ll*w`|| {. .};  
      ubound=uu*w`|| {. .};
      bound=lbound//ubound;    
      blc=bound//constrain;
      call nlpqn(rc, v, "dist", w, optn, blc);
      return (rc);
   finish fitwithbound;

   success=1;
   maxIterReached=0; 
   multiplier=0.5;
   ll=&lower; uu=&upper;
   initialL=&lower; initialU=&upper; finalrc=0;
   absolutebound=1;
   /*  initial */      
   rc=fitwithbound(ll, uu); 
   lastuu=uu;
   lastrc=rc;

   *print "before searching upper bound";
   *print initialL initialU ll uu lastuu rc;  
   if (rc<=0) then 
   do;
     uu=(1+multiplier)*lastuu; 
     absolutebound=absolutebound<>lastuu;
   end;

   if (rc>0) then
     uu=multiplier*lastuu<>absolutebound; 

   *print "going to search upper bound";
   *print absolutebound uu lastuu;  

   continue=1;
   iter=0;
   do while (continue=1);      
     rc=fitwithbound(ll, uu);
     if (rc<=0) then absolutebound=absolutebound<>uu;
      temp=uu;
      max=uu<>lastuu;
      min=uu><lastuu;
      diff=max-min;   
      if (rc<0 & lastrc<0) then
        uu=(1+multiplier)*max;
      /* *(max-absolutebound)*multiplier;  */
      if (rc>0 & lastrc>0) then
        uu=min-multiplier*(min-absolutebound); 
      if (rc>0 & lastrc<0) then
        uu=min+(1-multiplier)*diff;
      if (rc<0 & lastrc>0) then 
        uu=min+multiplier*diff; 
      /* always go above absolutebound */
      uu=uu<>absolutebound;
      
      lastuu=temp;
      lastrc=rc;

      if (diff<&eps & rc>0) then 
        do;    
          uu=lastuu;
          continue=0;
        end;
      if (diff=0 & rc<0) then do;
         ll=lastll;
         continue=0;
         success=0;
       end;
      iter=iter+1;
      *print iter rc lastrc ll uu lastuu min max absolutebound diff;
      if (iter>&maxiter & continue=1) then do;
        continue=0;
        maxIterReached=1;
        success=0;
      end;
   end; /* end of do while (continue=1); for looking for upper bound */

   printiterl=iter;
   /* now search for a lower bound */
   if (success=1) then do;
     **** search for lower bound ****; 
     lastll=ll;
     lastrc=rc;
     absolutebound=1;
     ll=lastll+multiplier*(absolutebound-lastll);  /* lastll should be 0.01 and rc>0 */
     continue=1;
     
     *print "before searching lower bound";
     *print ll uu iter rc; 
     do while (continue=1);      
       rc=fitwithbound(ll, uu);
       if (rc<=0) then absolutebound=absolutebound><ll;
       *print rc; 
       temp=ll;
       max=ll<>lastll;
       min=ll><lastll;
       diff=max-min; 
       if (rc<0 & lastrc<0) then 
         ll=&lowestbound+multiplier*(min-&lowestbound); 
       if (rc>0 & lastrc>0) then 
         ll=max+multiplier*(absolutebound-max); 
       if (rc>0 & lastrc<0) then 
         ll=min+multiplier*diff; 
       if (rc<0 & lastrc>0) then 
         ll=min+(1-multiplier)*diff;
       /* always go under absolutebound */
       ll=ll><absolutebound;
       lastll=temp;
       lastrc=rc;
       *print iter rc lastrc ll uu lastll min max absolutebound diff;
       if (diff<&eps & rc>0) then do;
         ll=lastll;
         continue=0;
       end;
       if (diff=0 & rc<0) then do;
         ll=lastll;
         continue=0;
         success=0;
       end;
       iter=iter+1;
       if (iter>&maxiter+printiterl & continue=1) then do;
         continue=0;
         maxIterReached=1;
         success=0;
       end;
     end; /* end while (continue=1); searching for lower bound */
   end; /* end of search lower end */
   if (success=1) then do;
     *print "The lower bound=" ll "Upper bound=" uu;    
     lbound=ll*w`|| {. .};  
     ubound=uu*w`|| {. .};
     bound=lbound//ubound;    
     blc=bound//constrain;
     call nlpqn(rc, v, "dist", w0, optn, blc);
     Calw=v`;  
     create _tmpnewwt_ from CalW[colname="&calwt"]; 
     append from CalW;
   end;
   call symputx('lower', ll);
   call symputx('upper', uu);
   call symputx('printiter', iter);
   call symputx('maxIterReached', maxIterReached);
   call symputx('success', success);
%goto CALREPWTS_LOGIT;  
/* end of LU_LOGIT method */

/* create replicate weights */
%CALREPWTS_LINEAR:
   %if (&success=0) %then %goto EXIT;
   %if ((&negtiveCalW=1) & (&method=NONE))
      %then %goto CALREPWTS_LOGIT;
   data &out; merge _temp_&data _tmpnewwt_ ;   
   data &out; merge _obs_&data &out; by _obs_;
   %if (&calRepWt=0) %then %goto EXIT;
   %do i=1 %to &nReps; 
       data _temp_; set &out;
          keep _obs_ &controlvar &repwtPrefix&i ;
          if (&repwtPrefix&i <=0) then delete; 
          run;
          proc iml;
             use _temp_;
             read all var {&repwtPrefix&i} into w;
             read all var {&controlVar} into x;
             close _temp_; 
             T={&ctrltotal}`; 
             wx=w#x; 
             beta=ginv(X`*wx)*(T-X`*w);
             CalW=w+w#(X*beta);
             create _tmpnewwt_ from CalW[colname="&repwtPrefix&i"];
             append from CalW;
          quit;
       data _temp_; merge _temp_ _tmpnewwt_;
       data &out; merge &out _temp_; by _obs_; 
   %end;
   %goto EXIT;  

%CALREPWTS_TRUNLINEAR:
   %if (&success=0) %then %goto EXIT;
   data &out; merge _temp_&data _tmpnewwt_ ;   
   data &out; merge _obs_&data &out; by _obs_;
   %if (&calRepWt=0) %then %goto EXIT;
   %do i=1 %to &nReps; 
       %let thisRepGood=1;
       data _temp_; set &out;
          keep _obs_ &controlvar &repwtPrefix&i ;
          if (&repwtPrefix&i <=0) then delete; 
          run;
          proc iml;
            sucessRep=1;
            use _temp_;
            read all var {&repwtPrefix&i} into w;
            read all var {&controlVar} into x;
            close _temp_;
            T={&ctrltotal}`;   
            j=J(nrow(T), 1, 0) ;
            constrain=x`||j|| T;
            w0 = w;
            optn = {0 0};
            blc  = constrain;
            start dist(v) global(w);
               y = (v/w`-J(1, nrow(w), 1))`;
               f=sum(w#(y##2));
               return(f);
            finish dist;
            start fitwithbound(ll, uu) global(w, optn, constrain);
               lbound=ll*w`|| {. .};  
               ubound=uu*w`|| {. .};
               bound=lbound//ubound;    
               blc=bound//constrain;
               call nlpqn(rc, v, "dist", w, optn, blc);
               return (rc);
            finish fitwithbound;
            ll=&lower; uu=&upper;
            lbound=ll*w`|| {. .};  
            ubound=uu*w`|| {. .};
            bound=lbound//ubound;    
            blc=bound//constrain;
            call nlpqn(rc, v, "dist", w, optn, blc);
            *print rc;
            iter=0;
            if (rc>0) then do;
               Calw=v`;  
               create _tmpnewwt_ from CalW[colname="&repwtPrefix&i"]; 
               append from CalW;
               end; 
            else do;
               multiplier=0.5;
               continue=1;
               do while (continue=1);
                  if (ll<=&lowestbound) then uu=(1+multiplier)*uu;
                  else ll=&lowestbound+multiplier*(ll-&lowestbound);
                  rc=fitwithbound(ll, uu);
                  if (rc<=0) then do;
                     if (ll<=&lowestbound) then uu=(1+multiplier)*uu;
                     else ll=&lowestbound+multiplier*(ll-&lowestbound);
                     end;
                  else continue=0;
                  iter=iter+1;
                  if (iter>&maxiter & continue=1) then do;
                     continue=0;
                     sucessRep=0;
                  end; 
               end;
               if (sucessRep=1) then do;
                  lbound=ll*w`|| {. .};  
                  ubound=uu*w`|| {. .};
                  bound=lbound//ubound;    
                  blc=bound//constrain;
                  call nlpqn(rc, v, "dist", w, optn, blc);
                  if (rc>0) then do;
                     Calw=v`;  
                     create _tmpnewwt_ from CalW[colname="&repwtPrefix&i"]; 
                     append from CalW;
                  end;
                  else sucessRep=0;
               end;  
            end;
            *print "for &i rep" ll uu sucessRep rc iter;
            call symputx('thisRepGood', sucessRep);
            quit;
            %if (&thisRepGood=1) %then %do; 
               data _temp_; merge _temp_ _tmpnewwt_;
               data &out; merge &out _temp_; by _obs_; 
               %end; 
            %else %do;
               %let badrep = %eval(&badrep. + 1);
               data &out; set &out; drop &repwtPrefix&i; 
            %end;
   %end;    
   %goto EXIT; 

%CALREPWTS_LOGIT:
   %if (&success=0) %then %goto EXIT;
   data &out; merge _temp_&data _tmpnewwt_ ;   
   data &out; merge _obs_&data &out; by _obs_;
   %if (&calRepWt=0) %then %goto EXIT;
   %do i=1 %to &nReps;
       %let thisRepGood=1;
       data _temp_; set &out;
          keep _obs_ &controlvar &repwtPrefix&i ;
          if (&repwtPrefix&i <=0) then delete; 
          run;
          proc iml;
             sucessRep=1;
             use _temp_;
             read all var {&repwtPrefix&i} into w;
             read all var {&controlVar} into x;
             close _temp_;
             T={&ctrltotal}`; 
             j=J(nrow(T), 1, 0) ;
             constrain=x`||j|| T;
             w0 = w;
             optn = {0 0};
             blc  = constrain;
             start dist(v) global(w);
                y = (J(1, nrow(w), 1))+ (v/w`)##((log(v/w`)-J(1, nrow(w), 1)));
                f=sum(y);
                return(f);
             finish dist;
             start fitwithbound(ll, uu) global(w, optn, constrain);
                lbound=ll*w`|| {. .};  
                ubound=uu*w`|| {. .};
                bound=lbound//ubound;    
                blc=bound//constrain;
                call nlpqn(rc, v, "dist", w, optn, blc);
                return (rc);
             finish fitwithbound;
             ll=&lower; uu=&upper;
             lbound=ll*w`|| {. .};  
             ubound=uu*w`|| {. .};
             bound=lbound//ubound; 
             blc=bound//constrain;
             call nlpqn(rc, v, "dist", w0, optn, blc);
             iter=0;
             if (rc>0) then do;
               Calw=v`;  
               create _tmpnewwt_ from CalW[colname="&repwtPrefix&i"]; 
               append from CalW;
               end; 
            else do;
               multiplier=0.5;
               continue=1;
               do while (continue=1);
                  if (ll<=&lowestbound) then uu=(1+multiplier)*uu;
                  else ll=&lowestbound+multiplier*(ll-&lowestbound);
                  rc=fitwithbound(ll, uu);
                  if (rc<=0) then do;
                     if (ll<=&lowestbound) then uu=(1+multiplier)*uu;
                     else ll=&lowestbound+multiplier*(ll-&lowestbound);
                     end;
                  else continue=0;
                  iter=iter+1;
                  if (iter>&maxiter & continue=1) then do;
                     continue=0;
                     sucessRep=0;
                  end; 
               end;
               if (sucessRep=1) then do;
                  lbound=ll*w`|| {. .};  
                  ubound=uu*w`|| {. .};
                  bound=lbound//ubound;    
                  blc=bound//constrain;
                  call nlpqn(rc, v, "dist", w, optn, blc);
                  if (rc>0) then do;
                     Calw=v`;  
                     create _tmpnewwt_ from CalW[colname="&repwtPrefix&i"]; 
                     append from CalW;
                  end;
                  else sucessRep=0;
               end;  
            end;
            *print "for &i rep" ll uu sucessRep rc iter;
            call symputx('thisRepGood', sucessRep);
            quit;
            %if (&thisRepGood=1) %then %do; 
               data _temp_; merge _temp_ _tmpnewwt_;
               data &out; merge &out _temp_; by _obs_; 
               %end; 
            %else %do;
               %let badrep = %eval(&badrep. + 1);
               data &out; set &out; drop &repwtPrefix&i; 
            %end;
   %end;    
   %goto EXIT; 

%EXIT: 

OPTIONS NOTES errors=max; 
ods listing;
 
%if ((&method=NONE) & (&success=1)) %then
  %do;
    %if (&switchToExp=0) %then 
        %put NOTE: The calibration weights &calwt are created by using the LINEAR method.;
    %else 
        %put WARNING: The calibration weights &calwt are created by using the EXPONENTIAL method, because the LINEAR method results negative calibration weights.; 
  %end;

%if ((&method=LINEAR) & (&negtiveCalW=1)) %then
  %put WARNING: The LINEAR method results negative calibration weights, consider use other methods or do not specify the method= parameter.;

%if ((&method=EXPONENTIAL)& (&success=1)) %then
  %put NOTE: The calibration weights &calwt are created by using the EXPONENTIAL method.;

%if ((&maxIterReached=1) & (&success=0)) %then %do;
   %if ((&method=TRUNLINEAR) & (&upperBSearchNeeded=1)
         & (&lowerBSearchNeeded=0)) %then   
   %put NOTE: The upper bound cannot be obtained after &maxiter iterations with specified LOWER=&lower.;
   %if ((&method=TRUNLINEAR) & (&upperBSearchNeeded=0)
         & (&lowerBSearchNeeded=1)) %then   
   %put NOTE: The lower bound cannot be obtained after &maxiter iterations with specified UPPER=&upper.;
   %if ((&method=TRUNLINEAR) & (&upperBSearchNeeded=1)
         & (&lowerBSearchNeeded=1)) %then   
   %put NOTE: The lower and upper bound cannot be obtained after &maxiter iterations.;
   %end;
   
%if (&success=0) %then %do;
    %if (((&method=TRUNLINEAR)|(&method=LOGIT)) & (&upperBSearchNeeded=0)
         & (&lowerBSearchNeeded=0)) %then       
      %put ERROR: Calibration weights cannot be constructed with the LOWER=&lower and UPPER=&upper bounds for the &method method.;
    %if (((&method=TRUNLINEAR)|(&method=LOGIT)) & (&upperBSearchNeeded=0)
         & (&lowerBSearchNeeded=1)) %then       
      %put ERROR: Calibration weights cannot be constructed with the UPPER=&upper for the &method method.;
    %if (((&method=TRUNLINEAR)|(&method=LOGIT)) & (&upperBSearchNeeded=1)
         & (&lowerBSearchNeeded=0)) %then       
      %put ERROR: Calibration weights cannot be constructed with the LOWER=&lower for the &method method.;    
      %PUT ERROR: Failed to construct calibration weights.;
    %end;
    
%if (&success=1) %then %do;
    %if (((&method=TRUNLINEAR)|(&method=LOGIT)) & (&upperBSearchNeeded=1)
          & (&lowerBSearchNeeded=0)) %then
        %put NOTE: After &printiter iterations, the upper bound is set to UPPER=&upper for the &method method.;
    %if (((&method=TRUNLINEAR)|(&method=LOGIT)) & (&upperBSearchNeeded=0)
          & (&lowerBSearchNeeded=1)) %then
        %put NOTE: After &printiter iterations, the lower bound is set to LOWER=&lower for the &method method.;
    %if (((&method=TRUNLINEAR)|(&method=LOGIT)) & (&upperBSearchNeeded=1)
          & (&lowerBSearchNeeded=1)) %then
        %put NOTE: After &printiter iterations, the lower bound is set to LOWER=&lower and the upper bound is set to UPPER=&upper for the &method method;
    %if ((&method=TRUNLINEAR)|(&method=LOGIT)) %then
        %put NOTE: The calibration weights &calwt are created by using the &method method with LOWER=&lower and UPPER=&upper bounds.;
    %if ((&calRepWt=1) & (&badrep>0)) %then 
        %put WARNING: &badrep replicates are removed from data set &out because the calibration process over corresponding replicate weights failed.;
%end;

    ods graphics on;
    ods listing;
%mend SurveyCalibrate;

