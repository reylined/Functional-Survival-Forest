# Functional-Survival-Forest
Implementation and simulation study of FSF

functions.r:  
Support functions for data simulation and `MFPCA` method.
Simulates longitudinal and survival data from multivariate joint model.

mFACEs_sim.r:  
Pipeline using `mFACEs` to extract features of longitudinal covariates
and RSF/Cox to fit the survival model

mFPCA_sim.r:  
Pipeline using `MFPCA` to extract features of longitudinal covariates
and RSF/Cox to fit the survival model

summary_AUC_BS.r:  
Calculates AUC and BS using output of mFACEs_sim.r and mFPCA_sim.r
