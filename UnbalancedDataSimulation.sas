* a sas program to randomly introduce missing data into previously made
simulation data sets;

*suppress the log file with the following statement;
options nonotes nosource nosource2 errors=0 nocenter;

*or permit the log file to print with the following statements to help debug program;
*options notes source source2 errors=20;

*suppress printing the output with the following statement;
ods listing close;

*or permit printing the output with the following statement to help debug program;
*ods listing;

*use the following option to print log output of logic statements and iteration info
during macro execution;
*options mprint symbolgen mlogic;

%let niter=1000;
%macro missing;
%do i = 1 %to &niter;

data e10g400it&i; 
infile "C:\Documents and Settings\Administrator\My Docume
nts\Correlations\GDDHtPt1DATA\e10g400it&i"; 
input env rep geno trait $ pheno;

data e4g250it&i; set e10g400it&i; if env < 5;if geno < 251;
data e2g75it&i; set e4g250it&i; if env < 3; if geno < 76;

*sample 5%, 15%, or 25% missing data points...
distribute the missing data three different ways:
0 type = missing whole plots (both traits) always.
50 type = 50% missing plots, 50% missing data randomly
distributed among remaining data points.
100 type = 100% missing data points randomly distributed among
data points;

*first do 0 type by making a multivariate data set from original data set;
data g4; set e4g250it&i; if trait = "GDD"; gdd = pheno;
data y4; set e4g250it&i; if trait = "HT"; ht = pheno;
data multiv4; merge g4 y4;

*add a random number to each observation to use in randomly sampling missing data points;
%let seed = 10000*&i;
rand = ranuni(&seed);
proc sort; by rand;
proc rank out = multiv4r; var rand; ranks rank;

*introduce z5, z15, and z30 as the new phenotypic values corresponding with
5%, 15%, and 30% missing data;
data e405a; set multiv4r; if rank > 100;
data e405g; set e405a; z05 = gdd; trait = "GDD"; drop gdd ht;
data e405y; set e405a; z05 = ht; trait = "HT"; drop gdd ht;
data e405; set e405g e405y; proc sort; by env rep geno trait;

data e4015a; set multiv4r; if rank > 300;
data e4015g; set e4015a; z015 = gdd; trait = "GDD"; drop gdd ht;
data e4015y; set e4015a; z015 = ht; trait = "HT"; drop gdd ht;
data e4015; set e4015g e4015y; proc sort; by env rep geno trait;

data e4025a; set multiv4r; if rank > 500;
data e4025g; set e4025a; z025 = gdd; trait = "GDD"; drop gdd ht;
data e4025y; set e4025a; z025 = ht; trait = "HT"; drop gdd ht;
data e4025; set e4025g e4025y; proc sort; by env rep geno trait;

**********;

data g2; set e2g75it&i; if trait = "GDD"; gdd = pheno;
data y2; set e2g75it&i; if trait = "HT"; ht = pheno;
data multiv2; merge g2 y2;

*add a random number to each observation to use in randomly sampling missing data points;
rand = ranuni(&seed+1);
proc sort; by rand;
proc rank out = multiv2r; var rand; ranks rank;

*introduce z5, z15, and z30 as the new phenotypic values corresponding with
5%, 15%, and 30% missing data;
data e205a; set multiv2r; if rank > 15;
data e205g; set e205a; z05 = gdd; trait = "GDD"; drop gdd ht;
data e205y; set e205a; z05 = ht; trait = "HT"; drop gdd ht;
data e205; set e205g e205y; proc sort; by env rep geno trait;

data e2015a; set multiv2r; if rank > 45;
data e2015g; set e2015a; z015 = gdd; trait = "GDD"; drop gdd ht;
data e2015y; set e2015a; z015 = ht; trait = "HT"; drop gdd ht;
data e2015; set e2015g e2015y; proc sort; by env rep geno trait;

data e2025a; set multiv2r; if rank > 75;
data e2025g; set e2025a; z025 = gdd; trait = "GDD"; drop gdd ht;
data e2025y; set e2025a; z025 = ht; trait = "HT"; drop gdd ht;
data e2025; set e2025g e2025y; proc sort; by env rep geno trait;

***************************************************************************;

*do "100 type" by randomly sampling missing data from among the 2n observations
where n = number of plots;
data e4100; set e4g250it&i; 

*add a random number to each observation to use in randomly sampling missing data points;
rand = ranuni(&seed+1);
proc sort; by rand;
proc rank out = e4100r; var rand; ranks rank;

data e41005; set e4100r; if rank > 200;z1005 = pheno; drop pheno;
proc sort; by env rep geno trait;
data e410015; set e4100r; if rank > 600;z10015 = pheno; drop pheno;
proc sort; by env rep geno trait;
data e410025; set e4100r; if rank > 1000;z10025 = pheno; drop pheno;
proc sort; by env rep geno trait;

**********;

data e2100; set e2g75it&i; 

*add a random number to each observation to use in randomly sampling missing data points;
rand = ranuni(&seed+1);
proc sort; by rand;
proc rank out = e2100r; var rand; ranks rank;

data e21005; set e2100r; if rank > 30;z1005 = pheno; drop pheno;
proc sort; by env rep geno trait;
data e210015; set e2100r; if rank > 90;z10015 = pheno; drop pheno;
proc sort; by env rep geno trait;
data e210025; set e2100r; if rank > 150;z10025 = pheno; drop pheno;
proc sort; by env rep geno trait;

***************************************************************************;

*do "50 type" by first sampling half of the missing data points from plots
from the multivariate data sets, after dropping those plots, recreate the tall
type data set and randomly sample the remaining data points completely at random;

data e4505a; set multiv4r; if rank > 50;
data e4505g; set e4505a; z505 = gdd; trait = "GDD"; drop gdd ht;
data e4505y; set e4505a; z505 = ht; trait = "HT"; drop gdd ht;
data e45015a; set multiv4r; if rank > 150;
data e45015g; set e45015a; z5015 = gdd; trait = "GDD"; drop gdd ht;
data e45015y; set e45015a; z5015 = ht; trait = "HT"; drop gdd ht;
data e45025a; set multiv4r; if rank > 250;
data e45025g; set e45025a; z5025 = gdd; trait = "GDD"; drop gdd ht;
data e45025y; set e45025a; z5025 = ht; trait = "HT"; drop gdd ht;

data e4505x; set e4505g e4505y; 
rand = ranuni(&seed+1);
proc sort; by rand;
proc rank out = e4505xr; var rand; ranks rank;
data e4505; set e4505xr; if rank > 100;
proc sort; by env rep geno trait;

data e45015x; set e45015g e45015y; 
rand = ranuni(&seed+1);
proc sort; by rand;
proc rank out = e45015xr; var rand; ranks rank;
data e45015; set e45015xr; if rank > 300;
proc sort; by env rep geno trait;

data e45025x; set e45025g e45025y; 
rand = ranuni(&seed+1);
proc sort; by rand;
proc rank out = e45025xr; var rand; ranks rank;
data e45025; set e45025xr; if rank > 500;
proc sort; by env rep geno trait;
**********;

data e2505a; set multiv2r; if rank > 7;
data e2505g; set e2505a; z505 = gdd; trait = "GDD"; drop gdd ht;
data e2505y; set e2505a; z505 = ht; trait = "HT"; drop gdd ht;
data e25015a; set multiv2r; if rank > 22;
data e25015g; set e25015a; z5015 = gdd; trait = "GDD"; drop gdd ht;
data e25015y; set e25015a; z5015 = ht; trait = "HT"; drop gdd ht;
data e25025a; set multiv2r; if rank > 37;
data e25025g; set e25025a; z5025 = gdd; trait = "GDD"; drop gdd ht;
data e25025y; set e25025a; z5025 = ht; trait = "HT"; drop gdd ht;

data e2505x; set e2505g e2505y; 
rand = ranuni(&seed+1);
proc sort; by rand;
proc rank out = e2505xr; var rand; ranks rank;
data e2505; set e2505xr; if rank > 16;
proc sort; by env rep geno trait;

data e25015x; set e25015g e25015y; 
rand = ranuni(&seed+1);
proc sort; by rand;
proc rank out = e25015xr; var rand; ranks rank;
data e25015; set e25015xr; if rank > 46;
proc sort; by env rep geno trait;

data e25025x; set e25025g e25025y; 
rand = ranuni(&seed+1);
proc sort; by rand;
proc rank out = e25025xr; var rand; ranks rank;
data e25025; set e25025xr; if rank > 76;
proc sort; by env rep geno trait;

***************************************************************************;

*now merge together all the data sets, the final data set for each parameter
combination should have 3 x 3 = 9 new variables;

data e4final; merge e405 e4015 e4025 e4505 e45015 e45025 e41005 e410015 e410025;
by env rep geno trait;
data e2final; merge e205 e2015 e2025 e2505 e25015 e25025 e21005 e210015 e210025;
by env rep geno trait;

*store data sets on files on hard drive;

data e4final; set e4final;
file "C:\Documents and Settings\Administrator\My Docume
nts\Correlations\GDDHtPt2DATA\e4g250miss\e4g250missit&i"; 
put env rep geno trait z05 z015 z025 z505 z5015 z5025 z1005 z10015 z10025;
data e2final; set e2final;
file "C:\Documents and Settings\Administrator\My Docume
nts\Correlations\GDDHtPt2DATA\e2g75miss\e2g75missit&i"; 
put env rep geno trait z05 z015 z025 z505 z5015 z5025 z1005 z10015 z10025;

*end iterative do-loop;
%end;
%mend missing;

%missing;
run;