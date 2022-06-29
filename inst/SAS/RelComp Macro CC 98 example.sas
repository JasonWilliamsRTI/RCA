/*
This version of the Reliabile Components Analysis program recreates tables from Cliff and Caruso (1998).

/*

%RELCOMP(dataset, ,# of components retained, reliability cutoff value for retained components, variables, reliabilities)

If the cutoff value is set to 1 then the number of components value is used, otherwise retention is based on the reliability cutoff.

*/


/***********************************************************************************************/
/***********************************************************************************************/
/*SECTION 1. Data simulation for data corresponding to that reported in Cliff and Caruso (1998)

Data are generated using the correlations in Table x, page X as population values.

/***********************************************************************************************/
/***********************************************************************************************/

proc iml;

/****************************************************/
/* 1) correlations from C and C to use as pop values
/****************************************************/

CORR=
{
1	0.5	0.83	0.66	0.68	0.62	0.55	0.55	0.6	0.45	0.43,
0.5	1	0.54	0.55	0.48	0.4	0.41	0.36	0.46	0.29	0.36,
0.83	0.54	1	0.65	0.74	0.68	0.57	0.58	0.58	0.46	0.45,
0.66	0.55	0.65	1	0.59	0.55	0.51	0.46	0.6	0.43	0.4,
0.68	0.48	0.74	0.59	1	0.61	0.56	0.52	0.49	0.4	0.45,
0.62	0.4	0.68	0.55	0.61	1	0.58	0.49	0.55	0.45	0.44,
0.55	0.41	0.57	0.51	0.56	0.58	1	0.55	0.6	0.63	0.43,
0.55	0.36	0.58	0.46	0.52	0.49	0.55	1	0.55	0.45	0.39,
0.6	0.46	0.58	0.6	0.49	0.55	0.6	0.55	1	0.66	0.42,
0.45	0.29	0.46	0.43	0.4	0.45	0.63	0.45	0.66	1	0.25,
0.43	0.36	0.45	0.4	0.45	0.44	0.43	0.39	0.42	0.25	1
}
;

/*************************/
/* 2) name, means, and N 
/*************************/
varNames = {information digitspan vocab arith comp simil piccomp picarr block objass digsym};
Mean = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; /* population means */
N = 1000000; /* sample size */


/*************************/
/* 3) simulate data
/*************************/

call randseed(123);
X = RandNormal(N, Mean, Corr); /* x is a 1000 x 3 matrix */
SampleMean = mean(X);
SampleCorr = corr(X);

/*************************/
/* 4) verify data
/*************************/
print corr samplecorr[colname=varNames rowname=VarNames];
print samplemean     [colname=varNames];

/*****************************/
/* 5) output dataset for RCA
/*****************************/
create CC from X[colname=varNames]; append from X; close;
quit;

proc contents data=CC; run;
proc corr data=CC; run;


/*******************************************************************************************/
/*******************************************************************************************/
/*END Data simulation
/*******************************************************************************************/
/*******************************************************************************************/


/**************************************************************************************************************************/
/**************************************************************************************************************************/
/**************************************************************************************************************************/
/**************************************************************************************************************************/
/**************************************************************************************************************************/
/**************************************************************************************************************************/


/*******************************************************************************************/
/*******************************************************************************************/
/*Reliaible components estimation
/*******************************************************************************************/
/*******************************************************************************************/

/*macro */

%macro RelComp(dataset, retain, criterion, varlist, reliab);

/********************************/
/*create dataset of correlations*/
/********************************/

ODS EXCLUDE ALL;
proc corr data=&dataset; var &varlist;
ods output PearsonCorr=corrs;
run;

data corrs; set corrs; keep &varlist; run;
ODS EXCLUDE NONE;



/***********************************************************************************************/
/***********************************************************************************************/
/*IML part 1 
/***********************************************************************************************/
/***********************************************************************************************/

PROC IML;

/***********************************************************************************/
/*1) read correlation matrix into IML*/
/***********************************************************************************/
use corrs;
read all 
into CORR;

/*************************/
/* 2) reliability vector */
/*************************/
RELROW = 
{
&reliab
}
;
REL=RELROW`;

/*************************/
/* 3) name vector */
/*************************/
NAMEROW = 
{
&varlist
}
;
VARNAMES=NAMEROW`;


/***********************************************************************************/
/* 4) Create RSTAR, matrix with corrs on off-diagonal and reliabilities on diagonal*/
/***********************************************************************************/

DIAG=diag(CORR);
REL2=DIAG(REL);
DIFF=DIAG-REL2;
RSTAR = CORR-DIFF;

/******************************************************************/
/* 5) generalized eigenvalues - equal to component reliabilities
/*    generalized eigenvectors - vectors of weights that are not
/*    normed to have sum of squared elements = 1
/******************************************************************/
call geneig(RELIABILITES, _WEIGHTS_NONORM, RSTAR, CORR); 


/**********************************************************************************************************/
/*6) normalize weights so that sum of squared elements in each column = 1
/*   NOTE: this is only needed to get comparable scale weights as in Table 3 of C&C p 298
/**********************************************************************************************************/

/*number of items placed into a macro variable*/
varcount=nrow(_WEIGHTS_NONORM);
call symputx("nvar", varcount);

/*norm weight vectors*/
_WEIGHTS_NORMED = _WEIGHTS_NONORM;

%do i = 1 %to &NVAR;
_WEIGHTS_NORMED[ ,&i] = _WEIGHTS_NORMED[ ,&i]/(NORM(_WEIGHTS_NONORM[ , &i]));
%end;

/***************************************************************************/
/*7) subset weights matrices to retained compoents 
     get loadings by multiplying (de normed) weights by correlation matrix
/***************************************************************************/
/*create vector of cutoff values*/
CRIT = j(&NVAR, 1, &CRITERION);
CRITCOMP = RELIABILITES >= CRIT;
critrel=sum(CRITCOMP);
call symputx("critrel", critrel);

%let NCOMPS = &RETAIN; 
%if &critrel gt 0 %then %let NCOMPS = &critrel;

_WEIGHTS_NONORMSUB = _WEIGHTS_NONORM[,1:&NCOMPS];
_WEIGHTS_NORMEDSUB = _WEIGHTS_NORMED[,1:&NCOMPS];

LOADSSUB = CORR*_WEIGHTS_NONORMSUB;

/***************************************************************************/
/*8) write out datasets for later processing                               */
/***************************************************************************/

create _varnames 			from VARNAMES; append from VARNAMES; close;
create _reliabilities 		from RELIABILITES; append from RELIABILITES; close;
create _corr 				from CORR; append from CORR; close;
create _rstar 				from RSTAR; append from RSTAR; close;
create _weights_nonorm 		from _WEIGHTS_NONORM; append from _WEIGHTS_NONORM; close;
create _weights_normed 		from _WEIGHTS_NORMED; append from _WEIGHTS_NORMED; close;
create _weights_normedsub 	from _WEIGHTS_NORMEDSUB; append from _WEIGHTS_NORMEDSUB; close;
create _weights_nonormsub 	from _WEIGHTS_NONORMSUB; append from _WEIGHTS_NONORMSUB; close;
create _loadssub 			from LOADSSUB; append from LOADSSUB; close;

quit;


/***********************************************************************************************/
/***********************************************************************************************/
/*IML part 1 ends
/***********************************************************************************************/
/***********************************************************************************************/

/*********************************************/
/*9) PROC FACTOR for varimax rotation matrix
/*********************************************/

/*rotate loadings */
proc transpose data=_loadssub out=_PATTERN; *id varname; run;
data _PATTERN; set _PATTERN end=last; run;
data _PATTERN2(type=factor); set _PATTERN; _type_ ='PATTERN'; run;

ODS EXCLUDE ALL;
proc factor data=_PATTERN2 rotate=varimax; ods output OrthTrans=_vartrans OrthRotFactPat=_ROTPAT; run;
proc print data=_vartrans; run;
ODS EXCLUDE NONE;

/***********************************************************************************************/
/***********************************************************************************************/
/*IML part 2 starts
/***********************************************************************************************/
/***********************************************************************************************/

PROC IML;

/*********************************************/
/*10) read back in necessary datasets
/*********************************************/
use _corr; read all into CORR;
use _rstar; read all into RSTAR;
use _weights_nonorm; 		read all into _WEIGHTS_NONORM;    
use _weights_normed; 		read all into _WEIGHTS_NORMED;    
use _weights_normedsub; 	read all into _WEIGHTS_NORMEDSUB; 
use _weights_nonormsub; 	read all into _WEIGHTS_NONORMSUB; 
use _loadssub; 				read all into LOADSSUB; 
use _vartrans; 				read all into VARTRANS; 

_LOADSROT = LOADSSUB*VARTRANS;
_WEIGHTS_NORMEDSUB_ROT = _WEIGHTS_NORMEDSUB*VARTRANS;
_WEIGHTS_NONORMSUB_ROT = _WEIGHTS_NONORMSUB*VARTRANS; /*july 15 - this should the matrix for WGT = WG*T in C and C*/

/***********************************************************/
/*11) Normalize rotated weights, comparable to Table 5 C&C
/***********************************************************/

%let NCOMPS = &RETAIN; 
%if &critrel gt 0 %then %let NCOMPS = &critrel;

WGTROTNORMED= _WEIGHTS_NORMEDSUB_ROT;

%do k = 1 %to &NCOMPS;
WGTROTNORMED[ ,&k] = WGTROTNORMED[ ,&k]/(NORM(_WEIGHTS_NORMEDSUB_ROT[ , &k]));
%end;

create wgtrotnonorm		from _WEIGHTS_NONORMSUB_ROT; append from _WEIGHTS_NONORMSUB_ROT; close;
create wgtrotnormed		from WGTROTNORMED; append from WGTROTNORMED; close;
create _loadsrot 		from _LOADSROT; append from _LOADSROT; close;

/***************************************************************************/
/*12) component reliabilites from rotated weights                              */
/***************************************************************************/

RELROT = WGTROTNORMED`*RSTAR*WGTROTNORMED*inv(WGTROTNORMED`*CORR*WGTROTNORMED);
RELROTVEC = vecdiag(RELROT);
create _reliabilitiesrot	from RELROTVEC; append from RELROTVEC; close;

quit;

/***********************************************************************************************/
/***********************************************************************************************/
/*IML part 2 ends
/***********************************************************************************************/
/***********************************************************************************************/


/*****************************************************************/
/* 13) Generate output for reliabilities table and weights table  */
/*****************************************************************/

/*******************/
/*table of weights */
/*******************/

data _varnames; set _varnames; rename COL1 = Varname; run;

/*************/
/*normed     */
/*************/
data _weights_normed; merge _varnames _weights_normed; nob = _N_; 
	%do i=1 %to &nvar ;
	rename    COL&i. = COMP&i.WGT;
	%end;
run;

/*create macro variables for original variables names*/
%do i = 1 %to &nvar;
proc sql noprint;
   select varname
      	into :var&i 
      	from _weights_normed
		where nob=&i;
   quit;
run;
%end;

/*************/
/*non normed */
/*************/
data _weights_nonorm; merge _varnames _weights_nonorm; nob = _N_; 
	%do i=1 %to &nvar ;
	rename    COL&i. = COMP&i.WGT;
	%end;
run;

/*create macro variables for original variables names*/
%do i = 1 %to &nvar;
proc sql noprint;
   select varname
      	into :var&i 
      	from _weights_nonorm
		where nob=&i;
   quit;
run;
%end;

/*****************************************************/
/*rotated non normed weights subset to k components  */
/*june 28                                            */
/****************************************************/

data wgtrotnonorm; merge _varnames wgtrotnonorm; nob = _N_; 
	%do i=1 %to &NCOMPS ;
	rename    COL&i. = COMP&i.ROTNONORM_WGT;
	%end;
run;

/************************************************/
/*rotated normed weights subset to k components */
/************************************************/

data wgtrotnormed; merge _varnames wgtrotnormed; nob = _N_; 
	%do i=1 %to &NCOMPS ;
	rename    COL&i. = COMP&i.ROT_WGT;
	%end;
run;

/*create macro variables for original variables names*/
%do i = 1 %to &nvar;
proc sql noprint;
   select varname
      	into :var&i 
      	from wgtrotnormed
		where nob=&i;
   quit;
run;
%end;

/******************************************/
/*rotated loadings subset to k components */
/******************************************/

data _loadsrot; merge _varnames _loadsrot; nob = _N_; 
	%do i=1 %to &NCOMPS ;
	rename    COL&i. = COMP&i.ROT_LOAD;
	%end;
run;

/*create macro variables for original variables names*/
%do i = 1 %to &nvar;
proc sql noprint;
   select varname
      	into :var&i 
      	from _loadsrot
		where nob=&i;
   quit;
run;
%end;

/*******************/
/*reliabilities     */
/*******************/
data _reliabilitiesrot; set _reliabilitiesrot; rename COL1 = RCA_rel_rot; run;
data _reliabilities; merge _reliabilities _reliabilitiesrot; rename COL1 = RCA_reliability; 
compnum=_N_;
RelComponent=cats("Component",compnum);
drop compnum;
run;
data _reliabilities; format RelComponent RCA_reliability RCA_rel_rot; set _reliabilities;
Label 
RCA_reliability = 'Unrotated reliabilties of all components'
RCA_rel_rot = 'Reliabilties of retained rotated components'; 
run;

/*****************************/
/*minor clean up for printing*/
/*****************************/
data _weights_nonorm_out; set _weights_nonorm;
drop nob; rename varname = Variable; run;

data _weights_normed_out; set _weights_normed;
drop nob; rename varname = Variable; run;

data wgtrotnonorm_out; set wgtrotnonorm;
drop nob; rename varname = Variable; run;

data wgtrotnormed_out; set wgtrotnormed;
drop nob; rename varname = Variable; run;

data _loadsrot_out; set _loadsrot;
drop nob; rename varname = Variable; run;


/*send to output*/
title1 'RCA component reliability estimates';
proc print label noobs data=_reliabilities; run;
title1;
title1 'Table of non-normed RCA weights for each component';
proc print noobs data=_weights_nonorm_out; run;
title1;
title1 'Table of normalized RCA weights for each component';
proc print noobs data=_weights_normed_out; run;
title1;
title1 'Table of orthogonally rotated RCA weights for retained components';
proc print noobs data=wgtrotnonorm_out; run;
title1;
title1 'Table of orthogonally rotated and normed RCA weights for retained components';
proc print noobs data=wgtrotnormed_out; run;
title1;
title1 'Table of orthogonally rotated loadings for retained components';
proc print noobs data=_loadsrot_out; run;
title1;

/*********************************************************************/
/*14)create single row of weights for composite creation in main data*/
/*********************************************************************/
data _ROTPAT; 
*set wgtrotnormed_out; 
*set _loadsrot_out; 
set wgtrotnonorm_out;
rename Variable=Varname; run;

data allweights; set _NULL_; run;
	%do i = 1 %to &NCOMPS;
	data temp; set _ROTPAT; keep COMP&i.: Varname; run;
	proc transpose data=temp out=_part&i suffix=_RCAwgt&i;
	id varname;
	var COMP&i.:;
	run;
	data allweights; merge allweights _part&i; key = 1; drop _NAME_; run;
	%end;

data &dataset; set &dataset; key=1; run;
data allweights ; set allweights; key=1; run;
data &dataset.RCA; merge &dataset allweights; by key; run;

/******************************************************************/
/*15)store weights names as macro strings                          */
/******************************************************************/

%do j = 1 %to &NCOMPS;
	proc transpose data=_part&j out=_flipped&j; run;
	data _flipped&j; set _flipped&j; nob=_N_; run;

		%do i = 1 %to &nvar;
		proc sql noprint;
		   select _NAME_
		      	into :wgtcomp&i._&j
		      	from _flipped&j
				where nob=&i;
		   quit;
		%end;
run;
%end;

/******************************************************************/
/*16)standardize Xs for composites                                */
/******************************************************************/

proc stdize data=&dataset.RCA out=_temp method=STD SPREFIX=std OPREFIX ; var &varlist; run;

data &dataset.RCA ; set _temp; run;

/******************************************************************/
/*17)make composites                                               */
/******************************************************************/
/*proc standard data=&dataset.RCA mean=0 stddev=1 out=temp; var &varlist; run;
data &dataset.RCA ; set temp; run;
*/

data &dataset.RCA; set &dataset.RCA;
%do j = 1 %to &NCOMPS;
	%do i = 1 %to &nvar;
	prod&j._&i = std&&var&i * &&wgtcomp&i._&j; 
	%end;
RELCOMP&j = sum(of prod&j.: );
%end;
run;

title1 'Component Score means and correlations';
proc corr data= &dataset.RCA; var RELCOMP1-RELCOMP&NCOMPS; run;
title1 ;


/*cleanup*/
proc datasets library=work;
   delete _: _flipped: _part: allweights temp rca_weights varnames corrs WGTROTNONORM WGTROTNORMED 
;
run;

/******************************/
/* 17)The End
/******************************/



%mend;

%RelComp(CC,5, .70, information digitspan vocab arith comp simil piccomp picarr block objass digsym, .90 .89 .96 .86 .84 .82 .85 .76 .89 .68 .86);
