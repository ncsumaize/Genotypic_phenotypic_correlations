* a sas program to analyze data sets with different amounts and different distributions 
of missing data with both REML and MANOVA;

*suppress the log file with the following statement;
options nonotes nosource nosource2 errors=0 nocenter;

*or permit the log file to print with the following statements to help debug program;
*options notes source source2 errors=20 center;

*suppress printing the output with the following statement;
ods listing close;

*or permit printing the output with the following statement to help debug program;
*ods listing;

*use the following option to print log output of logic statements and iteration info
during macro execution;
*options mprint symbolgen mlogic;

*or suppress log output of logic statements and iteration info
during macro execution;
options nomprint nosymbolgen nomlogic;

*make empty data sets to store ouput;

data glme4; 
data glme2; 
data mxe4; 
data mxe2; 
 
*********************************************************************************************
*DEFINE MACRO MANOVA to use GLM MANOVA to obtain variance and covariance components estimates;
%macro manova(dataset, pheno);

*make a multivariate data set from the original "tall data set";
data g&dataset; set &dataset; if trait = "GDD"; gdd = &pheno;
data h&dataset; set &dataset; if trait = "HT"; HT = &pheno;
data multiv; merge g&dataset h&dataset; by env rep geno;

proc glm data = multiv outstat = glmout;
class env rep geno;
model gdd HT = env rep(env) geno geno*env;
random env rep(env) geno geno*env;
manova h = env rep(env) geno geno*env;
ods output ExpectedMeanSquares = EMS;
quit;
ods output close;

*transform expected mean squares output to all capital letters to 
ensure correct processing of coefficients;
data EMS; set EMS;
source = upcase(Source);
ExpectedMeanSquare = upcase(ExpectedMeanSquare);

*get the coefficients on expected mean squares;
data coefGE; set EMS;if upcase(Source) = "ENV*GENO";
lGE = length(ExpectedMeanSquare);
nGE = lGE - 27;
coefGEGE = substr(ExpectedMeanSquare, 14, nGE);
data coefG; set EMS; if upcase(Source) = "GENO";
coefGE1 = substr(ExpectedMeanSquare, 14, 6);
coefGEG = compress(coefGE1, "VAR(ENV*GENO)");
lcoefGEG = length(coefGEG);
coefG1 = substr(ExpectedMeanSquare, 31 + lcoefGEG, 6);
coefGG = compress(coefG1, "VAR(GENO)");
data coefs; merge coefGE coefG; keep coefGEGE coefGEG coefGG;

*get the mean squares and cross products;
data glmout; set glmout; if _type_ = "ERROR" or _type_ = "SS3";
_source_ = upcase(_source_);

data error1; set glmout; 
if _source_ = "ERROR" and UPCASE(_NAME_) = "GDD";
	CovE = HT/DF; VE1 = SS/DF; dfe = DF;
data error2 ; set glmout;
if _source_ = "ERROR" and UPCASE(_NAME_) = "HT" ;
	VE2 = SS/DF; 
data ge1; set glmout;
if _source_ = "ENV*GENO" and UPCASE(_NAME_) = "GDD";
	MCPGE = HT/DF; MSGE1 = SS/DF; dfge = DF;
data ge2; set glmout;
if _source_ = "ENV*GENO" and UPCASE(_NAME_) = "HT";
	MSGE2 = SS/DF; 
data g1; set glmout;
if _source_ = "GENO" and UPCASE(_NAME_) = "GDD";
	MCPG = HT/DF;  MSG1 = SS/DF; dfg = DF;
data g2; set glmout;
if _source_ = "GENO" and UPCASE(_NAME_) = "HT";
	MSG2 = SS/DF; 

data statsout; merge error1 error2 GE1 GE2 G1 G2 coefs;

*solve for variance and covariance components estimates using method of moments;
CovGE = (MCPGE - CovE)/coefGEGE;
VGE1 = (MSGE1 - VE1)/coefGEGE;
VGE2 = (MSGE2 - VE2)/coefGEGE;
CovG = (MCPG - (coefGEG*CovGE) - CovE)/coefGG;
VG1 = (MSG1 - (coefGEG*VGE1) - VE1)/coefGG;
VG2 = (MSG2 - (coefGEG*VGE2) - VE2)/coefGG;
CovP = CovG + CovGE + CovE;
VP1 = VG1 + VGE1 + VE1;
VP2 = VG2 + VGE2 + VE2;
rp = CovP/(sqrt(VP1*VP2));

*compute variances and covariances of variance and covariance estimates;
VarVE1 = 2*(VE1**2)/(dfe+2);
VarVE2 = 2*(VE2**2)/(dfe+2);
VarCovE = ((VE1*VE2)+(CovE**2))/(dfe+2);
CovVE1VE2 = 2*(CovE**2)/(dfe+2);
CovVE1CovE = 2*VE1*CovE/(dfe+2);
CovVE2CovE = 2*VE2*CovE/(dfe+2);
VarVGE1 = ((1/coefGEGE)**2)*((2*(MSGE1**2)/(dfge+2))+VarVE1);
VarVGE2 = ((1/coefGEGE)**2)*((2*(MSGE2**2)/(dfge+2))+VarVE2);
VarCovGE = ((1/coefGEGE)**2)*((((MSGE1*MSGE2)+(MCPGE**2))/(dfge+2))+ VarCovE);
CovVGE1VGE2 = ((1/coefGEGE)**2)*((2*(MCPGE**2)/(dfge+2))+CovVE1VE2);
CovVGE1CovGE = ((1/coefGEGE)**2)*(((2*MSGE1*MCPGE)/(dfge+2))+CovVE1CovE);
CovVGE2CovGE = ((1/coefGEGE)**2)*(((2*MSGE2*MCPGE)/(dfge+2))+CovVE2CovE);
CovVGE1VE1 = -(1/coefGEGE)*VarVE1;
CovVGE1VE2 = -(1/coefGEGE)*CovVE1VE2;
CovVGE1CovE = -(1/coefGEGE)*CovVE1CovE;
CovVGE2VE1 = CovVGE1VE2;
CovVGE2VE2 = -(1/coefGEGE)*VarVE2;
CovVGE2CovE = -(1/coefGEGE)*CovVE2CovE;
CovCovGEVE1 = CovVGE1CovE;
CovCovGEVE2 = CovVGE2CovE;
CovCovGECovE = -(1/coefGEGE)*VarCovE;
*define the ratio of the coeficient on GE varcomp in the E(MS) of genotype vs. in the 
E(MS) of G*E, name it cr;
cr = coefGEG/coefGEGE;
VarVG1 = ((1/coefGG)**2)*((2*(MSG1**2)/(dfg+2))+(2*(cr**2)*(MSGE1**2)/(dfge+2))+(((1-cr)**2)*VarVE1));
VarVG2 = ((1/coefGG)**2)*((2*(MSG2**2)/(dfg+2))+(2*(cr**2)*(MSGE2**2)/(dfge+2))+(((1-cr)**2)*VarVE2));
VarCovG = ((1/coefGG)**2)*((((MSG1*MSG2)+(MCPG**2))/(dfg+2))+(((cr**2)*((MSGE1*MSGE2)+(MCPGE**2)))/(dfge+2))+(((1-cr)**2)*VarCovE));
CovVG1Vg2 = ((1/coefGG)**2)*(((2*(MCPG**2))/(dfg+2))+((2*(cr**2)*(MCPGE**2))/(dfge+2))+(((1-cr)**2)*CovVE1VE2));
CovVG1CovG = ((1/coefGG)**2)*(((2*MSG1*MCPG)/(dfg+2))+((2*(cr**2)*(MSGE1*MCPGE))/(dfge+2))+(((1-cr)**2)*CovVE1CovE));
CovVG2CovG = ((1/coefGG)**2)*(((2*MSG2*MCPG)/(dfg+2))+((2*(cr**2)*(MSGE2*MCPGE))/(dfge+2))+(((1-cr)**2)*CovVE2CovE));
CovVG1VGE1 = (1/(coefGG*coefGEGE))*(((-2*cr*(MSGE1**2))/(dfge+2))+((1-cr)*VarVE1));
CovVG1VGE2 = (1/(coefGG*coefGEGE))*(((-2*cr*(MCPGE**2))/(dfge+2))+((1-cr)*CovVE1VE2));
CovVG1CovGE = (1/(coefGG*coefGEGE))*(((-2*cr*(MSGE1*MCPGE))/(dfge+2))+((1-cr)*CovVE1CovE));
CovVG2VGE1 = CovVG1VGE2;
CovVG2VGE2 = (1/(coefGG*coefGEGE))*(((-2*cr*(MSGE2**2))/(dfge+2))+((1-cr)*VarVE2));
CovVG2CovGE = (1/(coefGG*coefGEGE))*(((-2*cr*(MSGE2*MCPGE))/(dfge+2))+((1-cr)*CovVE2CovE));
CovCovGVGE1 = (1/(coefGG*coefGEGE))*(((-2*cr*MSGE1*MCPGE)/(dfge+2))+((1-cr)*CovVE1CovE));
CovCovGVGE2 = (1/(coefGG*coefGEGE))*(((-2*cr*MSGE2*MCPGE)/(dfge+2))+((1-cr)*CovVE2CovE));
CovCovGCovGE = (1/(coefGG*coefGEGE))*(((-cr*((MSGE1*MSGE2)+(MCPGE**2)))/(dfge+2))+((1-cr)*VarCovE));
CovVG1VE1 = (1/coefGG)*(-(1-cr)*VarVE1);
CovVG1VE2 = (1/coefGG)*(-(1-cr)*CovVE1VE2);
CovVG1CovE = (1/coefGG)*(-(1-cr)*CovVE1CovE);
CovVG2VE1 = CovVG1VE2;
CovVG2VE2 = (1/coefGG)*(-(1-cr)*VarVE2);
CovVG2CovE = (1/coefGG)*(-(1-cr)*CovVE2CovE);
CovCovGVE1 = CovVG1CovE;
CovCovGVE2 = CovVG2CovE;
CovCovGCovE = (1/coefGG)*(-(1-cr)*VarCovE);

*IF EITHER GENETIC VARIANCE IS ESTIMATED TO BE LESS THAN OR EQUAL TO ZERO, 
NO CORRELATION ESTIMATE IS POSSIBLE, SO FORCE ESTIMATE AND ITS SE TO BE ZERO...;

if VG1 le 0 then
call SYMPUT("ZERO", "VG1");
else if VG2 le 0 then
call SYMPUT("ZERO", "VG2");
else call SYMPUT("ZERO", "NONE");

data statsout; set statsout;
%if &ZERO = VG1 %then %do;
rg = 0;
*create a fake variable Vg1b = 1 that will be used to make SE calculations function, this will allow
program to continue running, but later in IML step we replace the SE by zero in this case;
Vg1b=1;
Vg2b = VG2;
VarVG1 = 0;
CovVG1Vg2 = 0;
CovVG1CovG = 0;
CovVG1VGE1 = 0;
CovVG1VGE2 = 0;
CovVG1CovGE = 0;
CovVG1VE1 = 0;
CovVG1VE2 = 0;
CovVG1CovE = 0;
%end;

%else %if &ZERO = VG2 %then %do;
rg = 0;
Vg1b = Vg1;
Vg2b=1;
VarVG2 = 0;
CovVG1Vg2 = 0;
CovVG2CovG = 0;
CovVG2VGE1 = 0;
CovVG2VGE2 = 0;
CovVG2CovGE = 0;
CovVG2VE1 = 0;
CovVG2VE2 = 0;
CovVG2CovE = 0;
%end;

%else %do;
rg = CovG/(sqrt(VG1*VG2));
Vg1b = Vg1;
Vg2b = VG2;
%end;

proc iml;
*take data from statsout and read into vectors to be used to create matrices for calculations,
note that each vector will have one row and number of columns equal to the
number of variables. the order of variables is determined by their order in the 
'read' statements in iml NOT by the order in the original data set!;
 
*MAKE A VECTOR WITH GENETIC CORRELATION ESTIMATE;
use statsout; read var{rg} into rg; 

*MAKE A VECTOR WITH PHENOTYPIC CORRELATION ESTIMATE;
use statsout; read var{rp} into rp;

*MAKE A VECTOR WITH GENETIC VARIANCE AND COVARIANCE ESTIMATES;
use statsout; read var{VG1b CovG VG2b} into dg1; 

*MAKE A D VECTOR FOR COMPUTING STANDARD ERROR OF RG;
dg2 = vecdiag(inv(diag(dg1)));
*add zeroes to dg1 vector;
zero6 = {0 0 0 0 0 0}`;
dg3 = dg2//zero6;
*multiply (elementwise) dg3 by constant vector m;
m = {-0.5 1 -0.5 -0.5 1 -0.5 -0.5 1 -0.5}`;
dg = m#dg3; 

*MAKE A VECTOR WITH PHENOTYPIC VARIANCE AND COVARIANCE ESTIMATES;
use statsout; read var{VP1 CovP VP2} into dp1; 

*MAKE A D VECTOR FOR COMPUTING STANDARD ERROR OF RP;
dp2 = dp1||dp1||dp1;
dp3 = vecdiag(inv(diag(dp2)));
*multiply (elementwise) dp3 by constant vector m;
dp = m#dp3; 

*MAKE A VARIANCE-COVARIANCE MATRIX OF ESTIMATORS, "C";
use statsout; read var{
VarVG1 CovVG1CovG CovVG1Vg2 
CovVG1VGE1 CovVG1CovGE CovVG1VGE2 
CovVG1VE1 CovVG1CovE CovVG1VE2}
 into C1; 
use statsout; read var{
CovVG1CovG VarCovG CovVG2CovG 
CovCovGVGE1 CovCovGCovGE CovCovGVGE2
CovCovGVE1 CovCovGCovE CovCovGVE2}
 into C2; 
use statsout; read var{
CovVG1VG2 CovVG2CovG VarVG2 
CovVG2VGE1 CovVG2CovGE CovVG2VGE2 
CovVG2VE1 CovVG2CovE CovVG2VE2}
 into C3; 
use statsout; read var{
CovVG1VGE1 CovCovGVGE1 CovVG2VGE1 
VarVGE1 CovVGE1CovGE CovVGE1VGE2 
CovVGE1VE1 CovVGE1CovE CovVGE1VE2}
into C4; 
use statsout; read var{
CovVG1CovGE CovCovGCovGE CovVG2CovGE 
CovVGE1CovGE VarCovGE CovVGE2CovGE 
CovCovGEVE1 CovCovGECovE CovCovGEVE2} 
into C5; 
use statsout; read var{
CovVG1VGE2 CovCovGVGE2 CovVG2VGE2 
CovVGE1VGE2 CovVGE2CovGE VarVGE2 
CovVGE2VE1 CovVGE2CovE CovVGE2VE2}
into C6; 
use statsout; read var{ 
CovVG1VE1 CovCovGVE1 CovVG2VE1 
CovVGE1VE1 CovCovGEVE1 CovVGE2VE1 
VarVE1 CovVE1CovE CovVE1VE2}
into C7; 
use statsout; read var{ 
CovVG1CovE CovCovGCovE CovVG2CovE 
CovVGE1CovE CovCovGECovE CovVGE2CovE 
CovVE1CovE VarCovE CovVE2CovE} 
into C8; 
use statsout; read var{ 
CovVG1VE2 CovCovGVE2 CovVG2VE2
CovVGE1VE2 CovCovGEVE2 CovVGE2VE2
CovVE1VE2 CovVE2CovE VarVE2}
 into C9; 
*CREATE THE FULL C MATRIX BY CONCATENATING THE 9 COMPONENT VECTORS;
C = C1//C2//C3//C4//C5//C6//C7//C8//C9; 

*CALCULATE STANDARD ERRORS OF RG AND RP;
%if &ZERO = VG1 or &ZERO = VG2  %then %do;
SErg = 0;
%end;
%else %do;
SErg = Sqrt((rg**2)*dg`*C*dg);
%end;
SErp = Sqrt((rp**2)*dp`*C*dp);

*print rg SErg rp SErp;

*OUTPUT MATRIX COMPUTATIONS TO A DATA FILE;
final = rg||serg||rp||serp;
create og&dataset from final [colname = {'rg' 'serg' 'rp' 'serp'}];
append from final;

quit;

*COMBINE VARIANCE ESTIMATES AND VARIANCES OF VARIANCE ESTIMATORS TO COMPARE GLM TO MIXED 
(THE VARIANCES of VARIANCES SEEM TO BE BIGGER FROM MIXED;

data varvars; set statsout; 
iter = &i; pheno = "&pheno"; 
keep iter pheno VG1 VG2 CovG VarVG1 VarVG2 VarCovG;
data og2&dataset; merge og&dataset varvars;

*combine output from current iteration with that of previous iterations;
data check9;
ds = "&dataset";
if ds = :"e4it" then call symput("dset", "e4");
else if ds = :"e2it" then call symput ("dset", "e2");
run;

%if &dset = e4 %then %do;
	data glme4; length pheno $ 6; set glme4 og2&dataset;
	%end;
%else %if &dset = e2 %then %do;
	data glme2; length pheno $ 6; set glme2 og2&dataset;
	%end;

%mend manova;

*********************************************************************************************;
*DEFINE MACRO REML to use Multivariate REML to obtain variance and covariance components estimates;
%macro MREML(dataset, pheno);
proc mixed asycov data = &dataset; 
class trait env rep geno;
model &pheno = env(trait) rep(env*trait);
random trait/subject = geno type = un;
random trait/subject = geno*env type = un;
repeated trait/ sub = rep*geno(env) type = un;
ods output covparms = estmat; ods output asycov = covmat;

*Check to see if either variance component estimate is zero,
in which case, rg will be set to zero with zero SE;
data check1; set estmat; if upcase(subject) = "GENO"; if CovParm = "UN(1,1)";
VG1 = estimate; keep subject VG1;
data check2; set estmat; if upcase(subject) = "GENO"; if CovParm = "UN(2,2)";
VG2 = estimate; keep subject VG2;
data check; merge check1 check2;
if VG1 le 1E-10 then
call SYMPUT("REMLZERO", "VG1");
else if VG2 le 1E-10 then
call SYMPUT("REMLZERO", "VG2");
else call SYMPUT("REMLZERO", "NO");

data check; set check;

proc iml;
use estmat; read all into e;
use covmat; read all into cov;
* Note that SAS introduces an extra first column into the covariance matrix which must be removed;
C = cov(|1:nrow(cov), 2:ncol(cov)|); 
* Obtain genotypic and phenotypic covariance and variance components;
CovG = e(|2,1|);
VG1 = e(|1,1|);
VG2 = e(|3,1|);
CovP = CovG + e(|5,1|) + e(|8,1|);
VP1 = VG1 + e(|4,1|) + e(|7,1|);
VP2 = VG2 + e(|6,1|) + e(|9,1|);


* Create a module called "correl" that will estimate genotypic and phenotypic correlations 
and their standard errors;
start correl(C, CovG, VG1, VG2, CovP, VP1, VP2, RG, RP, SERG, SERP);
%if &REMLZERO = VG1 or &REMLZERO = VG2
%then %do; RG = 0;%end;
%else %do; RG = CovG/sqrt(VG1*VG2);%end;
*Make the derivative vector for rg, note that the order of the rows and columns of the variance
covariance matrix is VG1, CovG, VG2, VGE1, CovGE, VGE2, VError1, CovError, VError2;
%if &REMLZERO = VG1
%then %do; dg = (1//(1/CovG)//(-1/(2*VG2))//0//0//0//0//0//0;%end;
%else %if &REMLZERO = VG2
%then %do; dg = (-1/(2*VG1))//(1/CovG)//1//0//0//0//0//0//0;%end;
%else %do; dg = (-1/(2*VG1))//(1/CovG)//(-1/(2*VG2))//0//0//0//0//0//0;%end;
varrg = (RG**2)*dg`*C*dg; serg = sqrt(varrg); 
RP = CovP/sqrt(VP1*VP2);
*Make the derivate vector for rp;
d1p = -1/(2*VP1);
d2p = 1/CovP;
d3p = -1/(2*VP2);
dp= d1p//d2p//d3p//d1p//d2p//d3p//d1p//d2p//d3p;
varrp = (RP**2)*dp`*C*dp; serp = sqrt(varrp); 
finish correl;
call correl(C, CovG, VG1, VG2, CovP, VP1, VP2, RG, RP, SERG, SERP);

*get variances of variance components;
VarVG1 = C(|1,1|);
VarCoVG = C(|2,2|);
VarVG2 = C(|3,3|);

*OUTPUT MATRIX COMPUTATIONS TO A DATA FILE;
finalm = rg||serg||rp||serp||Vg1||CovG||Vg2||VarVG1||VarCoVG||VarVG2;
create om&dataset from finalm 
[colname = {'rg' 'serg' 'rp' 'serp' 'Vg1' 'CovG' 'Vg2' 'VarVG1' 'VarCovG' 'VarVG2'}];
append from finalm;

quit;

data om&dataset; set om&dataset;iter = &i; pheno = "&pheno";

*combine output from current iteration with that of previous iterations;
data check99;
ds = "&dataset";
if ds = :"e4it" then call symput("dst", "e4");
else if ds = :"e2it" then call symput ("dst", "e2");
run;

%if &dst = e4 %then %do;
	data mxe4; length pheno $ 6; set mxe4 om&dataset;
	%end;
%else %if &dst = e2 %then %do;
	data mxe2; length pheno $ 6; set mxe2 om&dataset;
	%end;

%mend MREML;

********************************************************************************************

*make a macro to do the simulation 1000 times;

%let niter=1000;
%macro simul;
%do i = 1 %to &niter;

data e4it&i; 
infile "C:\Documents and Settings\Administrator\My Docume
nts\Correlations\GDDHtUnbalancedDATA\e4g250miss\e4g250missit&i"; 
input env rep geno trait $ z05 z015 z025 z505 z5015 z5025 z1005 z10015 
z10025;

data e2it&i; 
infile "C:\Documents and Settings\Administrator\My Docume
nts\Correlations\GDDHtUnbalancedDATA\e2g75miss\e2g75missit&i"; 
input env rep geno trait $ z05 z015 z025 z505 z5015 z5025 z1005 z10015 
z10025;
*Invoke MANOVA and REML macros to analyze each data set;

%manova(e4it&i, z05);
%manova(e4it&i, z015);
%manova(e4it&i, z025);
%manova(e4it&i, z505);
%manova(e4it&i, z5015);
%manova(e4it&i, z5025);
%manova(e4it&i, z1005);
%manova(e4it&i, z10015);
%manova(e4it&i, z10025);
%manova(e2it&i, z05);
%manova(e2it&i, z015);
%manova(e2it&i, z025);
%manova(e2it&i, z505);
%manova(e2it&i, z5015);
%manova(e2it&i, z5025);
%manova(e2it&i, z1005);
%manova(e2it&i, z10015);
%manova(e2it&i, z10025);

%mreml(e4it&i, z05);
%mreml(e4it&i, z015);
%mreml(e4it&i, z025);
%mreml(e4it&i, z505);
%mreml(e4it&i, z5015);
%mreml(e4it&i, z5025);
%mreml(e4it&i, z1005);
%mreml(e4it&i, z10015);
%mreml(e4it&i, z10025);
%mreml(e2it&i, z05);
%mreml(e2it&i, z015);
%mreml(e2it&i, z025);
%mreml(e2it&i, z505);
%mreml(e2it&i, z5015);
%mreml(e2it&i, z5025);
%mreml(e2it&i, z1005);
%mreml(e2it&i, z10015);
%mreml(e2it&i, z10025);

%end; /*end iterative do-loop*/

%mend simul;

%simul;

*after completing all analyses, partition results according to
level and distribution of missing data;

%macro partition(dataset);

%macro partition2(pheno);
data &dataset.&pheno; set &dataset; if upcase(pheno) = "&pheno";

*store file on hard drive;
file "C:\Documents and Settings\Administrator\My Docume
nts\Correlations\GDDHtUnbalancedOUT\&dataset.&pheno"; 
put iter pheno rg SErg rp SErp Vg1 CovG Vg2 VarVG1 VarCovG VarVG2;
%mend partition2;

%partition2(Z05);
%partition2(Z015);
%partition2(Z025);
%partition2(Z505);
%partition2(Z5015);
%partition2(Z5025);
%partition2(Z1005);
%partition2(Z10015);
%partition2(Z10025);

%mend partition;

%partition(glme4);
%partition(glme2);
%partition(mxe4);
%partition(mxe2);
run;
