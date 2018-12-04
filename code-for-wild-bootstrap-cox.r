############### Functions for EE-bootstrap 
############### Example with cif estimation 
###########################################################################

wildCoxBand <- function(formula,data,xx,B=500,multiplier=c(2,3,1),Gsq=0)
 ### If Gsq=1, the weights G^2dN are used throughout in the wild bootstrap variances.
 ### multipliers 1=normal+1, 2=exp, 3=poisson
 ### B number of bootstraps 
 ### xx times where to look at baseline 
{ ## {{{
  
  rs <- 0 ## dN_i based resampling  inside timereg
  nj <- n <- nrow(data)

  ud <- tryCatch(
        cox.aalen(formula,data=data,
         rate.sim=rs,max.clust=NULL,max.timepoint.sim=NULL,basesim=1,Nit=100),
 	 error=function(x) NULL) 
                            
  if (!is.null(ud$gamma) && !is.na(ud$gamma)) {
    
    cumx <-  Cpred(ud$cum,xx)[,2]
    vcumx <- Cpred(ud$var.cum,xx)[,2]
    sdcumx <- vcumx^.5
    rvcumx <- Cpred(ud$robvar.cum,xx)[,2]
    sdrvcumx <- rvcumx^.5
    gem <- c(ud$gamma,ud$var.gamma^.5,ud$robvar.gamma^.5)
    ### 
    percen(ud$sim.testBeq0,0.95) ##, these are simulated sup Delta(t)/sdrvcumx(t)
    

    ###rbetadb <- rcumsdb <- c()
    rcumso <- rbeta <- rcums <-  rbeta.dir <- rcums.dir <- vcumdbs <- vcumdbs.dir <- c()
    successes <- 0
    failures <- 0
    
    repeat{
      rdd <- switch(multiplier[1],
                    "1" = {rnorm(nj)+1},
                    "2" = {rexp(nj)},
                    "3" = {rpois(nj,1)}
      )
      

      ud.EE <-tryCatch(
         cox.aalen(formula,data=data,rate.sim=rs,detail=0,caseweight=rdd,beta=ud$gamma,robust=0),error=function(x) NULL)
      

      if(is.null(ud.EE$gamma) || is.na(ud.EE$gamma) || any(is.na(Cpred(ud.EE$cum,xx)[,2]))) {
        
        failures <- failures + 1 
        
      }else{
        # sometimes ud.EE$gamma may have a huge (but non-NA) value and rcumx has only NaN entries (except for the first)
        
        successes <- successes + 1
        
        # for EE
        rcumx   <- Cpred(ud.EE$cum,xx)[,2]
        rbeta   <- rbind(rbeta,c(ud.EE$gamma))
        rcums   <- cbind(rcums,rcumx)
        rcumso   <- cbind(rcumso,ud.EE$cum[,2])
        
        vcumdbx <- Cpred(ud.EE$var.cum,xx)[,2]
        vcumdbs   <- cbind(vcumdbs,vcumdbx)
        
        
        # for "direct" resampling
        # here we use the centered variables, i.e. expectation 0.
        nj.dir <- length(ud$B.iid)

        rdd.dir <- switch(multiplier[1],
                      "1" = {rnorm(nj.dir)},
                      "2" = {rexp(nj.dir)-1},
                      "3" = {rpois(nj.dir,1)-1}
        )

        
        rcumx.dir <- Cpred( cbind(ud$cum[,1], sapply(ud$B.iid, identity) %*% rdd.dir) ,xx)[,2]
        rbeta.dir <- rbind(rbeta.dir, t(ud$gamma.iid) %*% rdd.dir)
        rcums.dir   <- cbind(rcums.dir,rcumx.dir)
        
        vcumdbx.dir <- Cpred( cbind(ud$cum[,1], sapply(ud$B.iid, identity)^2 %*% rdd.dir^2) ,xx)[,2]
        vcumdbs.dir <- cbind(vcumdbs.dir,vcumdbx.dir)
      }
      
      if(successes == B) break
      #if(successes == 1000) break
      #if(failures == 3000){
        #stop("Too many failures")
      #} 
      
      if (failures %% 100 == 99){
        print("Number of failures so far:")
        print(failures)
      } 
      
    }

    
    rcums.recovery <- rcums
    rcums <- na.omit(rcums)
    rbeta <- na.omit(rbeta)
    
    sdcumdb <- apply(rcums,1,sd)
    sdcumdb.dir <- apply(rcums.dir,1,sd)
    
    sdcumdb[sdcumdb==0] <- 1
    sdcumx[sdcumx==0] <- 1
    sdcumdb.dir[sdcumdb.dir==0] <- 1
    
    # a leading "c" stands for central, a leading "z" stands for studentized
    # ... for HW and EP bands, we don't need a studentization (at least: it is hidden in the weight function)
    ccums <- rcums-cumx
    zcums   <- (rcums-cumx)/sdcumdb
    ccums.dir <- rcums.dir
    zcums.dir <- rcums.dir / sdcumdb.dir
  
    cumx.inv <- 1/cumx
    # in order to not divide by 0
    cumx.inv[cumx.inv == 0] <- 1
    
    # In fact: here we use the same quantiles independent of log or not log. 
    # Therefore: A division by cumx is required in the definition of the log bands.
    pcumsdb.EE <- percen(apply(abs(zcums),2,max, na.rm=TRUE),0.95)
    # with "Hadamard"-derivative for the bootstrap version of the transformed NAE
    pcumsdb.EE.log <- percen(apply(abs(zcums),2,max, na.rm=TRUE),0.95)
    pcumsdb.dir <- percen(apply(abs(zcums.dir),2,max, na.rm=TRUE),0.95)
    # with "Hadamard"-derivative for the bootstrap version of the transformed NAE
    pcumsdb.dir.log <- percen(apply(abs(zcums.dir),2,max, na.rm=TRUE),0.95)
    

    band.EE <- cbind( cumx - sdcumx * pcumsdb.EE , cumx + sdcumx * pcumsdb.EE)
    band.EE.log <- cbind( cumx*exp(- pcumsdb.EE.log * sdcumx * cumx.inv)  ,cumx*exp( pcumsdb.EE.log * sdcumx * cumx.inv ))
    
    band.dir <- cbind(cumx - sdcumx * pcumsdb.dir , cumx + sdcumx * pcumsdb.dir)
    band.dir.log <- cbind( cumx*exp(- pcumsdb.dir.log * sdcumx * cumx.inv ) , cumx*exp( pcumsdb.dir.log * sdcumx * cumx.inv ))
    
    
    # Weight function for Hall-Wellner bands
    # for the estimator, and for their wild bootstrap versions
    HW <- 1 + n*sdcumx^2
    if(Gsq==0){
      HW.EE <- 1 + replicate(B, n*sdcumdb^2)
      HW.dir <- 1 + replicate(B, n*sdcumdb.dir^2)
    }else{
      # if Gsq==1, use the optional variation-type variance estimates
      HW.EE <- 1 + n*vcumdbs
      HW.EE <- pmax(HW.EE, 0.001)
      HW.dir <- 1 + n*vcumdbs.dir
      HW.dir <- pmax(HW.dir, 0.001)
    }
    
    HW.pcumsdb.EE <- percen(apply(abs(ccums / HW.EE),2,max, na.rm=TRUE),0.95)
    HW.pcumsdb.EE.log <- percen(apply(abs(ccums / HW.EE),2,max, na.rm=TRUE),0.95)
    HW.pcumsdb.dir <- percen(apply(abs(ccums.dir / HW.dir),2,max, na.rm=TRUE),0.95)
    HW.pcumsdb.dir.log <- percen(apply(abs(ccums.dir / HW.dir),2,max, na.rm=TRUE),0.95)
    band.HW.EE <- cbind( cumx - HW.pcumsdb.EE * HW , cumx + HW.pcumsdb.EE * HW)
    band.HW.EE.log <- cbind(cumx*exp(- HW.pcumsdb.EE.log * HW * cumx.inv ) , cumx*exp( HW.pcumsdb.EE.log * HW * cumx.inv ))
    band.HW.dir <- cbind(cumx - HW.pcumsdb.dir * HW ,cumx + HW.pcumsdb.dir * HW)
    band.HW.dir.log <- cbind(cumx*exp(- HW.pcumsdb.dir.log * HW * cumx.inv ), cumx*exp( HW.pcumsdb.dir.log * HW * cumx.inv ))
    
    
    # Weight function for equal precision bands
    # for the estimator, and for their wild bootstrap versions
    EP <- sqrt(n) * sdcumx
    if(Gsq==0){
      EP.EE <- replicate(B, sqrt(n)*sdcumdb)
      EP.dir <- replicate(B, sqrt(n)*sdcumdb.dir)
    }else{
      # if Gsq==1, use the optional variation-type variance estimates
      EP.EE <- sqrt(n * vcumdbs)
      EP.EE <- pmax(EP.EE, 0.001)
      EP.dir <- sqrt(n * vcumdbs.dir)
      EP.dir <- pmax(EP.dir, 0.001)
    }
    
    EP.pcumsdb.EE <- percen(apply(abs(ccums / EP.EE),2,max, na.rm=TRUE),0.95)
    EP.pcumsdb.EE.log <- percen(apply(abs(ccums / EP.EE),2,max, na.rm=TRUE),0.95)
    EP.pcumsdb.dir <- percen(apply(abs(ccums.dir / EP.dir),2,max, na.rm=TRUE),0.95)
    EP.pcumsdb.dir.log <- percen(apply(abs(ccums.dir / EP.dir),2,max, na.rm=TRUE),0.95)
    
    band.EP.EE <- cbind(cumx - EP.pcumsdb.EE * EP , cumx + EP.pcumsdb.EE * EP)
    band.EP.EE.log <- cbind(cumx*exp(- EP.pcumsdb.EE.log * EP * cumx.inv ) , cumx*exp( EP.pcumsdb.EE.log * EP * cumx.inv ))
    band.EP.dir <- cbind( cumx - EP * EP.pcumsdb.dir , cumx + EP * EP.pcumsdb.dir)
    band.EP.dir.log <-cbind( cumx*exp(- EP.pcumsdb.dir.log * EP * cumx.inv ) , cumx*exp( EP.pcumsdb.dir.log * EP * cumx.inv ))
  # band obtained via timereg (cox.aalen)
  band.timereg <- cbind(cumx-ud$conf.band*sdrvcumx ,cumx+ud$conf.band*sdrvcumx)
    
  } else {gem <- cumx <- vcumx <- rvcumx <- sdcumx <- sdrvcumx <- sdcumdb <- 
            cov.ud.band <- cov.EE <- cov.EE.log <- cov.dir <- cov.dir.log <- 
            HW.cov.ud.band <- HW.cov.EE <- HW.cov.EE.log <- HW.cov.dir <- HW.cov.dir.log <- 
            EP.cov.ud.band <- EP.cov.EE <- EP.cov.EE.log <- EP.cov.dir <- EP.cov.dir.log <- NULL;
    print("one iteration returned only NULLs")
  }
  
  ###gem <- c(gem,c(ud.EE$gamma,ud.EE$var.gamma^.5,ud.EE$robvar.gamma^.5,ud.EE$pval.Prop))
  return(list(x=xx,gem=gem,cum=cumx,vcum=vcumx,rvcum=rvcumx,
	      cumo=ud$cum[,2],betao=ud$gamma,
              rbeta=rbeta,rcums=rcums,rcumso=rcumso,
	      times=ud$cum[,1],  ### all the simuations within each
              sdcum=sdcumx, sdrvcumx=sdrvcumx,sdcumw=sdcumdb,
	      band.EE=band.EE, band.EE.log=band.EE.log,
              band.dir= band.dir, band.dir.log= band.dir.log,
	      band.HW.EE=band.HW.EE, band.HW.EE.log=band.HW.EE.log,
	      band.HW.dir=band.HW.dir, band.HW.dir.log=band.HW.dir.log,
	      band.EP.EE=band.EP.EE, band.EP.EE.log=band.EP.EE.log,
	      band.EP.dir=band.EP.dir, band.EP.dir.log=band.EP.dir.log,
	      band.timereg=band.timereg))
} ## }}} 

pred.cif.boot.timereg <- function(b1,b2,c1,c2,gplot=1)
{# {{{

times1 <- b1$times
times2 <- b2$times
coef1 <- b1$betao
coef2 <- b2$betao
###
bcums1 <- b1$rcumso
bcums2 <- b2$rcumso

where2 <- sindex.prodlim(times2,times1[-1],strict=TRUE)
cums2 <- c2$cum[,2]
cums2 <- cums2[where2]
cums1 <- c1$cum[,2]

bcums2 <- rbind(bcums2)[where2,]

n <- length(times1)
cif1 <- cumsum( exp(-cums1[-n]-cums2)*diff(cums1))

ccoef1 <- c1$gamma
ccoef2 <- c2$gamma
###

bcifs <- apply(exp(-bcums1[-n,]-bcums2)*apply(bcums1,2,diff),2,cumsum)

if (gplot==1) {
   matplot(times1[-1],bcifs,type="s",lwd=0.2)
   lines(times1[-1],cif1,type="s",lwd=2)
}

    cumx <- cif1
    ccums <- bcifs-cif1
    sdcumb <- apply(ccums,1,sd)
    zcums   <- ccums/sdcumb
###    ccums.dir <- rcums.dir
###    zcums.dir <- rcums.dir / sdcumdb.dir
  
    cumx.inv <- 1/cumx
    # in order to not divide by 0
    cumx.inv[cumx.inv == 0] <- 1
    
    # In fact: here we use the same quantiles independent of log or not log. 
    # Therefore: A division by cumx is required in the definition of the log bands.
    pcumsdb.EE <- percen(c(apply(abs(zcums),2,max, na.rm=TRUE)),0.95)
    # with "Hadamard"-derivative for the bootstrap version of the transformed NAE
    pcumsdb.EE.log <- percen(apply(abs(zcums),2,max, na.rm=TRUE),0.95)
###    pcumsdb.dir <- percen(apply(abs(zcums.dir),2,max, na.rm=TRUE),0.95)
    # with "Hadamard"-derivative for the bootstrap version of the transformed NAE
###    pcumsdb.dir.log <- percen(apply(abs(zcums.dir),2,max, na.rm=TRUE),0.95)
    

    band.EE <- cbind( cumx - sdcumb * pcumsdb.EE , cumx + sdcumb * pcumsdb.EE)
    band.EE.log <- cbind( cumx*exp(- pcumsdb.EE.log * sdcumb * cumx.inv)  ,cumx*exp( pcumsdb.EE.log * sdcumb * cumx.inv ))
###    cov.EE.log <- tt>= cumx*exp(- pcumsdb.EE.log * sdcumx * cumx.inv ) & tt<= cumx*exp( pcumsdb.EE.log * sdcumx * cumx.inv )
### 

    return(list(time=times1[-1],
	bcifs=bcifs,cif1=cif1,band.EE=band.EE,band.EE.log=band.EE.log))

}# }}}

library(mets)
###
data(TRACE)
data(tTRACE)
tTRACE$cage <- scale(tTRACE$age)
tTRACE$diabetes0 <- (tTRACE$diabetes==0)*1

ca <- cox.aalen(Surv(time,status!=0)~prop(diabetes)+prop(sex)+prop(cage),
	data=tTRACE,resample.iid=1,max.clust=NULL,max.timepoint.sim=NULL)
summary(ca)

llca <- wildCoxBand(Surv(time,status!=0)~prop(diabetes)+prop(sex)+prop(cage),tTRACE,seq(0,7,by=0.01),B=1000)
llca0 <- wildCoxBand(Surv(time,status!=0)~prop(diabetes0)+prop(sex)+prop(cage),tTRACE,seq(0,7,by=0.01),B=1000)

pca <- predict(ca,X=rbind(c(1),c(1)),Z=rbind(c(0,0,0),c(1,0,0)))
plot(pca,multiple=1,uniform=1,col=1:2,se=1,lty=1:2)
###
with(llca, matlines(x,exp(-band.EE.log),type="s",lty=3,col=3,lwd=2))
with(llca0, matlines(x,exp(-band.EE.log),type="s",lty=3,col=3,lwd=2))
### github version of mets 
### devtools::install_github("kkholst/mets")
### library(mets)
### mets:::plot.conf.region(llca$x,exp(-llca$band.EE.log),col=3)
### mets:::plot.conf.region(llca0$x,exp(-llca0$band.EE.log),col=3)

### restricted mean bootstrap 
coxrm <- restricted.residual.mean(ca,tau=5,x=rbind(c(0,0,0),c(1,0,0)),iid=1)
summary(coxrm)
plot(coxrm)
summary(coxrm)

#############################################################################
################### Cumulative incidence Wild Boostrap ######################
#############################################################################

data(bmt)
c1 <- cox.aalen(Surv(time,cause==1)~prop(platelet)+prop(age)+prop(tcell),bmt,robust=0)
c2 <- cox.aalen(Surv(time,cause==2)~prop(platelet)+prop(age)+prop(tcell),bmt,robust=0)

b1 <- bca1 <- wildCoxBand(Surv(time,cause==1)~prop(platelet)+prop(age)+prop(tcell),bmt,seq(0,100,by=1),B=1000)
b2 <- bca2 <- wildCoxBand(Surv(time,cause==2)~prop(platelet)+prop(age)+prop(tcell),bmt,seq(0,100,by=1),B=1000)

ud <- pred.cif.boot.timereg(b1,b2,c1,c2,gplot=1)
with(ud, plotl(time,cif1,ylim=c(0,1)))
matlines(ud$time,ud$band.EE,lty=1,col=2,type="s")
matlines(ud$time,ud$band.EE.log,lty=1,col=3,type="s")

