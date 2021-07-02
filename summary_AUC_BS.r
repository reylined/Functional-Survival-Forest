rm(list=ls())
setwd('PATH')

library(survival)
library(tdROC)
library(abind)
source('functions.R')

method = "mfpca"  #[mfpca, mfaces]
scenario = "none" #[none, interaction]
path = paste0("output/",method,"/",scenario,"/")

n.sim = 100      #number of simulation runs
n = 300          #sample size
n.train = 200    #n.test = n - n.train

#dynamic prediction information
obstime = seq(0,10,0.5)
Tstart = c(1,2,3,4) # landmark time for prediction
deltaT = c(1,2) # prediction windowns

row.names = cbind(rep(Tstart, each=length(deltaT)), rep(deltaT, length(Tstart)))

models = c("RSF","Cox")

listofTables = vector(mode="list", length=2)

i.model=1
for(model in models){

    #initalize empty matrix to store AUC/BS results
    auctrue.m = matrix(NA, nrow=nrow(row.names), ncol=n.sim)
    auc.m = matrix(NA, nrow=nrow(row.names), ncol=n.sim)
    BS.m = matrix(NA, nrow=nrow(row.names), ncol=n.sim)
    MSE.m = matrix(NA, nrow=nrow(row.names), ncol=n.sim) 
    
    for(i.run in 1:n.sim){
        load(paste0(path,method,"_",model,i.run,".rdata"))
        
        train.surv = surv[which(surv$id %in% train.id),]
        
        #estimate km curve for BS calculation
        km = survfit(Surv(time, event)~1, data=train.surv)
        survest = stepfun(km$time, c(1, km$surv))
        
        ith=0
        for(t in Tstart){
            for(dt in deltaT){
                ith = ith + 1
                tp = t + dt
                
                test.surv = surv[surv$id%in%DP.id[[ith]], ]
                N_vali = nrow(test.surv)

                # TRUE AUC
                roc = tdROC( X = 1-trueProb[[ith]], Y = timeEvent[[ith]]$time,
                             delta = timeEvent[[ith]]$event,
                             tau = tp, span = 0.05,
                             nboot = 0, alpha = 0.05,
                             n.grid = 1000, cut.off = 0.5)

                auctrue.m[ith, i.run] = roc$AUC$value
                
                # AUC
                roc = tdROC( X = 1-DP.prob[[ith]], Y = timeEvent[[ith]]$time, 
                             delta = timeEvent[[ith]]$event,
                             tau = tp, span = 0.05,
                             nboot = 0, alpha = 0.05,
                             n.grid = 1000, cut.off = 0.5)
                
                auc.m[ith, i.run] = roc$AUC$value
                
                # BRIER SCORE
                D = rep(0, N_vali)
                D[test.surv$time<=tp&test.surv$event==1] = 1
                pi = 1-DP.prob[[ith]]
                
                km_pts = survest(test.surv$time)/survest(t)
                W2 <- D/km_pts
                W1 <- as.numeric(test.surv$time>tp)/(survest(tp)/survest(t))
                W <- W1 + W2
                
                BS_pts <- W * (D - pi)^2
                BS.m[ith, i.run]  = sum(na.omit(BS_pts)) / N_vali
            }
        }
    }
    
    table = cbind(row.names,  apply(auctrue.m, 1, mean, na.rm=T),  apply(auc.m, 1, mean, na.rm=T),
                  apply(BS.m, 1, mean, na.rm=T), apply(MSE.m, 1, mean, na.rm=T))
    
    colnames(table) = c("t", "dt", "True.AUC", "AUC", "BS", "MSE")
    
    print(model)
    print(table)
    listofTables[[i.model]] = table
    i.model = i.model+1
}

table = cbind(listofTables[[1]][,1:5],listofTables[[2]][,4:5])
print(table)
save(table, file = paste0("output/output_",method,'_',scenario,'.rdata'))