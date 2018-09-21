* a sas program to draw random samples from a multivariate normal distribution
of two variables and estimate genetic and phenotypic correlations from each sample;

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

*make a macro to do the simulation 1000 times;

%let niter=1000;
%macro simul;
%do i = 1 %to &niter;

*generate pairs of random variables drawn from a multivariate normal distribution with 
known mean, variances, and covariances;
*for the full data set, there are 10 environments, 2 replications per environment, 
and 400 genotypes;

proc iml;
*first generate random environment effects with mean zero, variances of 7605 and 25, and covari-
ance of 448;
*note that this does not work because the specified var-cov matrix is not positive definite,
because the correlation between environment effects (448/sqrt(7605*25)) is > 1.  so,
scale the covariance to make the correlation 0.75, and the covariance should be 327;
seed=ceil(ranuni(0)*100000000);
mu0 = {0,0};
VCenv = {7605 327, 327 25};
call vnormal (env1, mu0, VCenv, 10, seed);
e = 1:10;
env = e`||env1;
create env from env [colname = {'Env' 'EnvX' 'EnvY'}];
append from env; 


*generate rep(env) effects with variances of 3 and 1 and covariance of -3;
*same problem as above, the matrix is not positive definite, because the correlation
between rep effects is < -1.  so, scale the covariance to -1.3 to make correlation -0.75;
VCrep = {3 -1.3,-1.3 1};
seed=seed+1;
call vnormal (rep1, mu0, VCrep, 20 , seed);
r = 1:20;
rep = r`||rep1;
create rep from rep [colname = {'Rep' 'RepX' 'RepY'}];
append from rep; 

*generate genotype effects with variances of 11001 and 11 and covariance of 116;

VCgeno = {11001 116, 116 11};
seed=seed+1;
*note - vnormal sometimes gives wrong output if output matrix is big,
SAS suggests calling it twice identically to solve problem...
i do NOT see this problem on Pentium 3, only Pentium 4 machines!;
call vnormal (geno1, mu0, VCgeno, 400 , seed);
call vnormal (geno1, mu0, VCgeno, 400 , seed);
g = 1:400;
geno = g`||geno1;
create geno from geno [colname = {'Geno' 'GenoX' 'GenoY'}];
append from geno; 

*generate GxE interaction effects with variances of 2319 and 13 and covariance 47;
VCge = {2319 47, 47 13};
seed=seed+1;
call vnormal (ge1, mu0, VCge, 4000, seed);
call vnormal (ge1, mu0, VCge, 4000, seed);
gbe = 1:4000;
gbe = gbe`||ge1;
create gbe from gbe [colname = {'GE' 'GEX' 'GEY'}];
append from gbe; 

*generate error effects with variances of 1452 and 31 with covariance of 20;
VCerr = {1452 20, 20 31}; 
seed=seed+1;
call vnormal (err, mu0, VCerr, 8000, seed);
call vnormal (err, mu0, VCerr, 8000, seed);
create err from err [colname = {'ErrX' 'ErrY'}];
append from err; 
quit;

*match up env, rep(env),  geno, genoxenv, err and overall effects;

data env; set env; 
data reps;
do env = 1 to 10;
	do rep = 1 to 2;output; 
	end; 
end;
data repli; merge rep reps;	
 
data gxe; 
do env = 1 to 10;
	do geno = 1 to 400;output; 
	end; 
end;
data gbye; merge gxe gbe;

data errs; 
do env = 1 to 10;
	do rep = 1 to 2;
		do geno = 1 to 400;output;
		end; 
	end; 
end;
data error; merge errs err; 

data envrep; merge env repli; by env;

data envreper; merge envrep error; by env rep;
proc sort; by geno;

data erg; merge envreper geno; by geno;
proc sort; by env geno;

data all; merge erg gbye; by env geno;
*add in overall mean effects of 1325 and 93;
GDD = 1325 + EnvX + RepX + GenoX + GEX + ErrX;
HT = 93 + EnvY + RepY + GenoY + GEY + ErrY;
*keep env rep geno GDD HT;

*now make "tall" data set for reml analysis;

data e10g400it&i; set all;
trait = "GDD"; pheno = gdd; output;
trait = "HT"; pheno = ht; output;
keep env rep geno trait pheno ;

*save data set as a text file on hard drive;
%macro storefile(dataset);
data &dataset; set &dataset;
file "C:\Documents and Settings\Administrator\My Docume
nts\Correlations\GDDHtBalancedDATA\&dataset"; 
*file "C:\Documents and Settings\Administrator\My Docume
nts\Jim\Correlations\&dataset"; 
put env rep geno trait pheno;
%mend storefile;

%storefile(e10g400it&i); 

%end; /*end iterative do-loop*/
%mend simul;
%simul;

run;