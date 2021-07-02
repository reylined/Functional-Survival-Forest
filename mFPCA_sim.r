rm(list=ls())
setwd('PATH')

library(MASS)
library(MFPCA)
library(randomForestSRC)
library(pec)
library(survival)
source('functions.r')

n.sim = 100     #number of simulation runs
n = 300         #sample size
n.train = 200   #n.test = n - n.train

#dynamic prediction information
obstime = seq(0,10,0.5) # longitudinal measurement time
Tstart = c(1,2,3,4) # landmark time for prediction
deltaT = c(1,2) # prediction windowns
argvals = obstime

scenario = "none"   #options: ["none","interaction"]


for(i.run in 1:n.sim){
    print(i.run)
    set.seed(123+i.run)
    
    # simulation
    sim.data = sim_mjm_linear(n, obstime=obstime, opt=scenario)
    long = sim.data$long  # longitudinal data
    surv = sim.data$surv  # survival data

    # split data to training and testing
    train.id = c(1:n.train)
    test.id = c((n.train+1):n)
    
    train.surv = surv[surv$id%in%c(1:n.train), ]
    test.surv = surv[!surv$id%in%c(1:n.train), ]
    
    # subject ids
    patID = surv$id
    nPat = length(patID)
    
    # transfer longitudinal outcomes from long to wide
    multivar = array(NA, c(n, length(obstime), 3))
    for(i in 1:nPat){
        visits = which(obstime %in% (long$obstime[long$id == patID[i]]))
        multivar[i,visits, 1] = long$Y1[long$id == patID[i]]
        multivar[i,visits, 2] = long$Y2[long$id == patID[i]]
        multivar[i,visits, 3] = long$Y3[long$id == patID[i]]
    }
    multivar.train = multivar[train.id, , ]
    
    
    # univariate FPCA via PACE
    Xi.train = L = phi.train = meanFun.train =  NULL
    x.tmp = NULL
    for(p in 1:3){
        tmp.ufpca = uPACE(multivar.train[,,p], argvals, nbasis=7)
        x.tmp = cbind(x.tmp,tmp.ufpca$scores[,1])
        Xi.train = cbind(Xi.train, tmp.ufpca$scores) # FPC scores
        L = c(L, dim(tmp.ufpca$scores)[2])
        phi.train[[p]] = t(tmp.ufpca$functions@X) # FPC eigenfunctions
        meanFun.train[[p]] = tmp.ufpca$mu@X # estimated mean functions
    }
    
    # multivariate FPCA
    mFPCA.train = mFPCA(Xi=Xi.train, phi=phi.train, p=3, L=L, I=n.train)
    rho.train = mFPCA.train$rho  #MFPC scores
    pve = mFPCA.train$pve
    psi = mFPCA.train$psi
    Cms = mFPCA.train$Cms
    
    colnames(rho.train) = paste0("rho.",(1:ncol(rho.train)))
    train.surv.temp = cbind(train.surv[,2:5],rho.train)
    
    # RSF Model
    rsf.fit = rfsrc(Surv(time,event)~., data=train.surv.temp,
                    ntree=1000, seed=i.run)
    
    # Cox Model
    cox.fit = coxph(Surv(time, event)~., data = train.surv.temp, model=T, x=T, y=T)

    
    # dynamic prediction
    DP.id = DP.prob.rsf = DP.prob.cox = DP.prob.truecox = timeEvent= trueProb = NULL
    ith = 0
    for(t in Tstart){
        tmp.id = test.surv[test.surv$time>t, "id"]  # subjects still event-free at landmark time 
        tmp.surv.data = test.surv[test.surv$time>t, ] # filter the data
        tmp.data = multivar[tmp.id, , ] # subset longitudinal outcomes
        tmp.data[,-c(1:which(t==obstime)),] = NA  # retain the measurements before landmark time
        
        # univariate FPC 
        Xi.test = NULL
        for(p in 1:3){
            tmp.ufpca = uPACE(multivar.train[,,p], argvals, tmp.data[,,p], nbasis = 7)
            Xi.test = cbind(Xi.test, tmp.ufpca$scores) # dynamic FPC scores for test subjects 
        }
        
        # estimate MFPC scores for test subjects
        rho.test = mfpca.score(Xi.test, Cms)
        colnames(rho.test) = paste0("rho.",(1:ncol(rho.test)))
        test.surv.temp = cbind(tmp.surv.data[,2:5],rho.test)
        
        # prediction for different time windowes
        for(dt in deltaT){
            ith = ith + 1
            DP.id[[ith]] = tmp.id
            timeEvent[[ith]] = tmp.surv.data[, c("time", "event")] # true event time and even indicator
            trueProb [[ith]] = tmp.surv.data$true.prob[, (which((t+dt)==obstime)-1)] # true risk 
            DP.prob.rsf[[ith]] = cond.prob.pec(rsf.fit, test.surv.temp, t, (t+dt))  # risk prediction
            DP.prob.cox[[ith]] = cond.prob.pec(cox.fit, test.surv.temp, t, (t+dt))  # risk prediction
        }
    }
    DP.prob = DP.prob.rsf
    save(sim.data, surv, train.id, trueProb, DP.id, DP.prob, timeEvent,
         file=paste(c("output/mfpca/",scenario,"/mfpca_rsf",i.run,".rdata"), collapse=""))

    DP.prob = DP.prob.cox
    save(sim.data, surv, train.id, trueProb, DP.id, DP.prob, timeEvent,
         file=paste(c("output/mfpca/",scenario,"/mfpca_cox",i.run,".rdata"), collapse=""))

}


