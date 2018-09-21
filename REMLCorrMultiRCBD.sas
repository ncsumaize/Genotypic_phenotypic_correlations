data one;
FILENAME DATA1 DDE "EXCEL|[DataForCorrelationPrograms.xls]Sheet1!R2C1:R1054C9";
INFILE DATA1 NOTAB DLM= '09'X DSD MISSOVER lrecl = 10240;
INPUT Env Rep Block Plot Ent Geno $ GDD Ht Yield;
proc print;
run;

*first, estimate variance components for each trait separately to compare to multivariate analysis below;
%macro varcomp(trait);
proc mixed data = one;
class env rep geno;
model &trait = ;
random env rep(env) geno geno*env;
*also check effect of setting environments and reps fixed on other variance components;
proc mixed data = one;
class env rep geno;
model &trait = env rep(env) ;
random geno geno*env;
run;
%mend;

%varcomp(GDD);
%varcomp(HT);
%varcomp(YIELD);

*restructure data set for multivariate reml analysis;
data two; length trait $ 5; set one;
trait = "GDD"; y = gdd; output;
trait = "HT"; y = ht; output;
trait = "YIELD"; y = yield; output;
drop gdd ht yield;

* analyze variables pair-wise;

%macro corr(trait1, trait2);
data traits; set two; if trait = "&trait1" or trait = "&trait2";
proc mixed asycov data = traits; 
class trait env rep geno;
model y = env(trait) rep(env*trait);
random trait/subject = geno type = un;
random trait/subject = geno*env type = un;
repeated trait/ sub = rep*geno(env) type = un;
ods output covparms = estmat; ods output asycov = covmat;
run;
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
RG = CovG/sqrt(VG1*VG2);
*Make the derivative vector for rg, note that the order of the rows and columns of the variance
covariance matrix is VG1, CovG, VG2, VGE1, CovGE, VGE2, VError1, CovError, VError2;
dg = (-1/(2*VG1))//(1/CovG)//(-1/(2*VG2))//0//0//0//0//0//0;
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
print "Genotypic Correlation Between &trait1 and &trait2";
print RG serg;
print "Phenotypic Correlation Between &trait1 and &trait2";
print RP serp;
quit; run;
%mend;

%corr(GDD,HT);
%corr(GDD,YIELD);
%corr(HT,YIELD);
run;
