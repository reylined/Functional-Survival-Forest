rm(list=ls())
setwd('PATH')

library(MASS)
library(mfaces)
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
argvals = seq(0,max(obstime),length=10*max(obstime)+1)

scenario = "none"   #options: ["none","interaction"]


for(i.run in 1:n.sim){
    print(i.run)
    set.seed(123+i.run)
    
    
    # simulation
    sim.data = sim_mjm_linear(n, obstime=obstime, opt=scenario)
    long = sim.data$long
    surv = sim.data$surv

    
    # split data to training and testing
    train.id = c(1:n.train)
    test.id = c((n.train+1):n)
    
    train.surv = surv[surv$id%in%train.id, ]
    test.surv = surv[surv$id%in%test.id, ]
    train.long = long[long$id%in%train.id, ]
    test.long = long[long$id%in%test.id, ]


    y1 = data.frame("subj"=train.long$id,"argvals"=train.long$obstime,"y"=train.long$Y1)
    y2 = data.frame("subj"=train.long$id,"argvals"=train.long$obstime,"y"=train.long$Y2)
    y3 = data.frame("subj"=train.long$id,"argvals"=train.long$obstime,"y"=train.long$Y3)
    
    y1 = y1[which(!is.na(y1$y)),]
    y2 = y2[which(!is.na(y2$y)),]
    y3 = y3[which(!is.na(y3$y)),]
    
    multivar.train = list(y1,y2,y3)
    
    # mFACEs
    mfaces.fit = mface.sparse(multivar.train, argvals.new = argvals, knots=7,
                        newdata = multivar.train, calculate.scores = TRUE)

    train.surv.tmp = cbind(train.surv[,2:5],mfaces.fit$scores$scores)
    
    # Random Survival Forest
    rsf.fit = rfsrc(Surv(time,event)~., data=train.surv.tmp,
                    ntree=1000, seed=i.run)
    
    # Cox Model
    cox.fit = coxph(Surv(time,event)~., data=train.surv.tmp, model=TRUE, x=TRUE, y=TRUE)
    
    # dynamic prediction
    DP.id = DP.prob = DP.prob.rsf = DP.prob.cox = timeEvent = trueProb = NULL
    ith = 0
    for(t in Tstart){
        tmp.surv.data = test.surv[test.surv$time>t, ] # subjects still event-free at landmark time
        tmp.id = tmp.surv.data[["id"]]
        test.long = test.long[which(test.long$id %in% tmp.id),]
        test.long = test.long[which(test.long$obstime<=t),]
        
        y1 = data.frame("subj"=test.long$id,"argvals"=test.long$obstime,"y"=test.long$Y1)
        y2 = data.frame("subj"=test.long$id,"argvals"=test.long$obstime,"y"=test.long$Y2)
        y3 = data.frame("subj"=test.long$id,"argvals"=test.long$obstime,"y"=test.long$Y3)
        
        y1 = y1[which(!is.na(y1$y)),]
        y2 = y2[which(!is.na(y2$y)),]
        y3 = y3[which(!is.na(y3$y)),]
        
        multivar.test = list("y1"=y1,"y2"=y2,"y3"=y3)
        
        mfaces.pred = predict(mfaces.fit, multivar.test, scores=TRUE)
        
        test.surv.tmp = cbind(tmp.surv.data[,2:5],mfaces.pred$scores$scores)
        
        # prediction for different time windowes
        for(dt in deltaT){
            ith = ith + 1
            DP.id[[ith]] = tmp.id
            timeEvent[[ith]] = tmp.surv.data[, c("time", "event")] # true event time and event indicator
            trueProb [[ith]] = tmp.surv.data$true.prob[, (which((t+dt)==obstime)-1)] # true risk 
            DP.prob.rsf[[ith]] = cond.prob.pec(rsf.fit, test.surv.tmp, t, (t+dt))  # risk prediction
            DP.prob.cox[[ith]] = cond.prob.pec(cox.fit, test.surv.tmp, t, (t+dt))  # risk prediction
        }
    }
    DP.prob = DP.prob.rsf
    save(sim.data, surv, train.id, trueProb, DP.id, DP.prob, timeEvent,
         file=paste(c("output/mfaces/",scenario,"/mfaces_rsf",i.run,".rdata"), collapse=""))
    
    DP.prob = DP.prob.cox
    save(sim.data, surv, train.id, trueProb, DP.id, DP.prob, timeEvent,
         file=paste(c("output/mfaces/",scenario,"/mfaces_cox",i.run,".rdata"), collapse=""))
}



