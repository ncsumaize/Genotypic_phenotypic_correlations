

EXAMPLE SAS PROGRAMS USED FOR SIMULATION STUDY DESCRIBED IN HOLLAND, J.B. "Estimating genotypic and phenotypic correlations and their standard errors using multivariate restricted maximum likelihood estimation in SAS Proc MIXED"

[SAS program used to generate full simulation data sets](FullSimulation.sas). 

(This program requires a folder with the following path: “C:\Documents and Settings\Administrator\My Documents\Correlations\GDDHtBalancedDATA” to be created.  This is the folder where the simulation data sets will be stored).

 

[SAS program used to analyze each balanced simulation data set with MANOVA (Proc GLM) and REML (Proc MIXED)](REMLvsMANOVABalanced.sas).

(This program requires two folders to be created:

"C:\Documents and Settings\Administrator\My Documents\Correlations\GDDHtUnbalancedDATA\e2g75miss"

and "C:\Documents and Settings\Administrator\My Documents\Correlations\GDDHtUnbalancedDATA\e4g250miss"

These two folders will store the unbalanced simulation data sets).

 

[SAS program used to make unbalanced data sets from full data sets](UnbalancedDataSimulation.sas). 

(This program requires a folder with the following path: "C:\Documents and Settings\Administrator\My Documents\Correlations\GDDHtBalancedOUT” to be created.  This is the folder where the parameter estimates from each analysis will be stored).

 

[SAS program used to analyze each unbalanced simulation data set with MANOVA (Proc GLM) and REML (Proc MIXED)](REMLvsMANOVAUnbalanced.sas).

(This program requires a folder with the following path: "C:\Documents and Settings\Administrator\My Documents\Correlations\GDDHtUnbalancedOUT” to be created.  This is the folder where the parameter estimates from each analysis will be stored).











