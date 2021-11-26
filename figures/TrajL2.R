#set.seed(14)
seed=12345678
setwd('R')
source('SimulationsFunctions.R')
source('Parameters.R')
source('Costs.R')

trunc=FALSE

###############################  Grille 1 #########################################""
distance="L2"
grille=642
		load(paste("RData/XGrid.RData",sep=""))
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )

step="all"
CVs=1; CVi=1;
	load(paste("RData/ProgDyn_",grille,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))				
################## Coût Simple


	set.seed(seed)
	xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
	Strat=list(dec=rep("non",N),vis=rep(0,N))
	yn=0
	Psin=as.matrix(c(1,rep(0,nOmega-2)))
	decision<-Decsimple[1,1]
	Strat$dec[1]<-Decisions[decision,1]
	Strat$vis[1]<-Decisions[decision,2]		
	RealTime=0
	n=1
	IndexGamma<-c(1)
	Dist=0

Modes<-c(1,0,0)
ModesEst<-c(1,0,0)

while(RealTime<H)
{
	if (n==1)
	{
		xnew=NewStep(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),j0="nonsense",Strat$dec[n],tjump=0,mjump=0,TsH=0)
	} 	else
	{
		xnew=NewStep(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n-1],Strat$dec[n],tjump,mjump,TsH)
	}	
	xn$traj[n+1]=xnew[2]
	xn$time[n+1]=xnew[3]
	xn$mode[n+1]=xnew[1]
	xn$realtime[n+1]=xnew[4]
	RealTime<-RealTime+as.numeric(Strat$vis[n])
	tjump=xnew[5]
	mjump=xnew[6]		
	TsH=xnew[7]


	if(xnew[1]==100)
	{
		ynew=0 ; Psin=c(Psin,rep(0,nOmega-1)); xn$realtime[n+2]=H; xn$traj[n+2]=D; Modes=rbind(Modes, c(0,0,0)); Modes=rbind(Modes, c(0,0,0)); RealTime<-H
	} else
	{	
		vmult=which(Time==Strat$vis[n])
		dmult=which(Treatment==Strat$dec[n])
		dind=(vmult-1)*3+dmult
		ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
		if (trunc)
			if (ynew<2) ynew=0			
		PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Psin[,n]
		eG<-evalGauss(ynew-xval)
		Num=PTemp*eG
		Den=eG%*%PTemp
		Psitemp=Num/as.vector(Den)
		Psin=cbind(Psin,Psitemp)
		NT=NNeighborT(Psitemp)[nOmega]
		indGamma=NT
		IndexGamma<-c(IndexGamma,indGamma)
		Mbelief<-Gamma[,indGamma]
		Dist=c(Dist,sqrt(sum((Psitemp-Gamma[,NT])^2)))
		Modes<-rbind(Modes,c(sum(Mbelief[1:nm0]),sum(Mbelief[(nm0+1):(nm0+nm1)]),sum(Mbelief[(nm0+nm1+1):(nm0+nm1+nm2)])))
		ModesEst<-rbind(ModesEst,c(sum(Psitemp[1:nm0]),sum(Psitemp[(nm0+1):(nm0+nm1)]),sum(Psitemp[(nm0+nm1+1):(nm0+nm1+nm2)])))
		indvisit=RealTime/delta
		if (RealTime<H)
		{
			decision<-Decsimple[indvisit+1,indGamma]
			Strat$dec[n+1]<-Decisions[decision,1]
			Strat$vis[n+1]<-Decisions[decision,2]
		} else 
		{
			Strat$dec[n+1]<-nd+1	
		}	
		yn=c(yn,ynew)
	}
	n<-n+1
}


	ll=which(xn$realtime==H)
	if (sum(xn$mode==100)>0) {ll=which(xn$mode==100)}
	Creal=round(RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"simple"),0)
	Cfilt=round(RCostFilt(IndexGamma,Strat$dec[1:ll],Strat$vis[1:ll],"simple"),0)



pdf(paste("figure/L2-",grille,".pdf",sep=""),width=15,height=10)
par(mfrow=c(5,1),mar=c(2,5,2,2))
plot(xn$realtime[1:ll],xn$traj[1:ll],type="b",col=ifelse(Strat$dec[1:ll]=="non","black",ifelse(Strat$dec[1:ll]=="a","green","red")),pch=ifelse(xn$mode[1:ll]==0,1,ifelse(xn$mode[1:ll]==1,2,3)),xlab="time (days)",ylab="real processus value",lwd=2,cex=2,cex.axis=2,cex.lab=2,xlim=c(0,H),ylim=c(0,10))
legend("top",legend=paste('Grid:',grille),bty="n",cex=2)
legend("topleft",legend=c("healthy","disease 1","disease 2"),pch=c(1:3),cex=2)
legend("topright",legend=c(paste("Real cost:",Creal),paste('Filtered cost:',Cfilt)),cex=2)
plot(xn$realtime[1:ll],yn[1:ll],type="b",col=ifelse(Strat$dec[1:ll]=="non","black",ifelse(Strat$dec[1:ll]=="a","green","red")),xlab="time (days)",ylab="observed data", lwd=2,cex=2,cex.axis=2,cex.lab=2,xlim=c(0,H),ylim=c(0,10))
legend("topleft",legend=c("no treatment","treatment a","treatment b"),col=c("black","green","red"),cex=2,lwd=2)
plot(xn$realtime[1:ll],ModesEst[1:ll,1],type="l",col="black",ylim=c(0,1),xlab="time (days)",ylab="estimated mode" ,lwd=2,cex=2,cex.axis=2,cex.lab=2,xlim=c(0,H))
legend("topleft",legend=c("mode 0","mode 1","mode 2"),col=c("black","green","red"),cex=2,lwd=2)
lines(xn$realtime[1:ll],ModesEst[1:ll,2],col="green",lwd=2,type="b")
lines(xn$realtime[1:ll],ModesEst[1:ll,3],col="red",lwd=2,type="b")
plot(xn$realtime[1:ll],Modes[1:ll,1],type="l",col="black",ylim=c(0,1),xlab="time (days)",ylab="projected mode" ,lwd=2,cex=2,cex.axis=2,cex.lab=2,xlim=c(0,H))
legend("topleft",legend=c("mode 0","mode 1","mode 2"),col=c("black","green","red"),cex=2,lwd=2)
lines(xn$realtime[1:ll],Modes[1:ll,2],col="green",lwd=2,type="b")
lines(xn$realtime[1:ll],Modes[1:ll,3],col="red",lwd=2,type="b")
plot(xn$realtime[1:ll],Dist,type="b",ylab="Dist to proj" ,lwd=2,cex.axis=2,cex.lab=2,xlim=c(0,H),ylim=c(0,0.5))
abline(h=0.35,col="grey")
dev.off()


IndexGamma

###############################  Grille 8 #########################################""
grille=2100
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )

step="all"
CVs=1; CVi=1;
	load(paste("RData/ProgDyn_",grille,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))				
################## Coût Simple


	set.seed(seed)
	xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
	Strat=list(dec=rep("non",N),vis=rep(0,N))
	yn=0
	Psin=as.matrix(c(1,rep(0,nOmega-2)))
	decision<-Decsimple[1,1]
	Strat$dec[1]<-Decisions[decision,1]
	Strat$vis[1]<-Decisions[decision,2]		
	RealTime=0
	n=1
	IndexGamma<-c(1)
	Dist=0

Modes<-c(1,0,0)
ModesEst<-c(1,0,0)

while(RealTime<H)
{
	if (n==1)
	{
		xnew=NewStep(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),j0="nonsense",Strat$dec[n],tjump=0,mjump=0,TsH=0)
	} 	else
	{
		xnew=NewStep(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n-1],Strat$dec[n],tjump,mjump,TsH)
	}	
	xn$traj[n+1]=xnew[2]
	xn$time[n+1]=xnew[3]
	xn$mode[n+1]=xnew[1]
	xn$realtime[n+1]=xnew[4]
	RealTime<-RealTime+as.numeric(Strat$vis[n])
	tjump=xnew[5]
	mjump=xnew[6]		
	TsH=xnew[7]


	if(xnew[1]==100)
	{
		ynew=0 ; Psin=c(Psin,rep(0,nOmega-1)); xn$realtime[n+2]=H; xn$traj[n+2]=D; Modes=rbind(Modes, c(0,0,0)); Modes=rbind(Modes, c(0,0,0)); RealTime<-H
	} else
	{	
		vmult=which(Time==Strat$vis[n])
		dmult=which(Treatment==Strat$dec[n])
		dind=(vmult-1)*3+dmult
		ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
		if (trunc)
			if (ynew<2) ynew=0			
		PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Psin[,n]
		eG<-evalGauss(ynew-xval)
		Num=PTemp*eG
		Den=eG%*%PTemp
		Psitemp=Num/as.vector(Den)
		NT=NNeighborT(Psitemp)[nOmega]
		Dist=c(Dist,sqrt(sum((Psitemp-Gamma[,NT])^2)))
		Psin<-cbind(Psin,Psitemp)
		indGamma=NT
		IndexGamma<-c(IndexGamma,indGamma)
		Mbelief<-Gamma[,indGamma]
		Modes<-rbind(Modes,c(sum(Mbelief[1:nm0]),sum(Mbelief[(nm0+1):(nm0+nm1)]),sum(Mbelief[(nm0+nm1+1):(nm0+nm1+nm2)])))
		ModesEst<-rbind(ModesEst,c(sum(Psitemp[1:nm0]),sum(Psitemp[(nm0+1):(nm0+nm1)]),sum(Psitemp[(nm0+nm1+1):(nm0+nm1+nm2)])))
		indvisit=RealTime/delta
		if (RealTime<H)
		{
			decision<-Decsimple[indvisit+1,indGamma]
			Strat$dec[n+1]<-Decisions[decision,1]
			Strat$vis[n+1]<-Decisions[decision,2]
		} else 
		{
			Strat$dec[n+1]<-nd+1	
		}	
		yn=c(yn,ynew)
	}
	n<-n+1
}


	ll=which(xn$realtime==H)
	if (sum(xn$mode==100)>0) {ll=which(xn$mode==100)}
	Creal=round(RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"simple"),0)
	Cfilt=round(RCostFilt(IndexGamma,Strat$dec[1:ll],Strat$vis[1:ll],"simple"),0)




pdf(paste("figure/L2-",grille,".pdf",sep=""),width=15,height=10)
par(mfrow=c(5,1),mar=c(2,5,2,2))
plot(xn$realtime[1:ll],xn$traj[1:ll],type="b",col=ifelse(Strat$dec[1:ll]=="non","black",ifelse(Strat$dec[1:ll]=="a","green","red")),pch=ifelse(xn$mode[1:ll]==0,1,ifelse(xn$mode[1:ll]==1,2,3)),xlab="time (days)",ylab="real processus value",lwd=2,cex=2,cex.axis=2,cex.lab=2,xlim=c(0,H),ylim=c(0,10))
legend("top",legend=paste('Grid:',grille),bty="n",cex=2)
legend("topleft",legend=c("healthy","disease 1","disease 2"),pch=c(1:3),cex=2)
legend("topright",legend=c(paste("Real cost:",Creal),paste('Filtered cost:',Cfilt)),cex=2)
plot(xn$realtime[1:ll],yn[1:ll],type="b",col=ifelse(Strat$dec[1:ll]=="non","black",ifelse(Strat$dec[1:ll]=="a","green","red")),xlab="time (days)",ylab="observed data", lwd=2,cex=2,cex.axis=2,cex.lab=2,xlim=c(0,H),ylim=c(0,10))
legend("topleft",legend=c("no treatment","treatment a","treatment b"),col=c("black","green","red"),cex=2,lwd=2)
plot(xn$realtime[1:ll],ModesEst[1:ll,1],type="l",col="black",ylim=c(0,1),xlab="time (days)",ylab="estimated mode" ,lwd=2,cex=2,cex.axis=2,cex.lab=2,xlim=c(0,H))
legend("topleft",legend=c("mode 0","mode 1","mode 2"),col=c("black","green","red"),cex=2,lwd=2)
lines(xn$realtime[1:ll],ModesEst[1:ll,2],col="green",lwd=2,type="b")
lines(xn$realtime[1:ll],ModesEst[1:ll,3],col="red",lwd=2,type="b")
plot(xn$realtime[1:ll],Modes[1:ll,1],type="l",col="black",ylim=c(0,1),xlab="time (days)",ylab="projected mode" ,lwd=2,cex=2,cex.axis=2,cex.lab=2,xlim=c(0,H))
legend("topleft",legend=c("mode 0","mode 1","mode 2"),col=c("black","green","red"),cex=2,lwd=2)
lines(xn$realtime[1:ll],Modes[1:ll,2],col="green",lwd=2,type="b")
lines(xn$realtime[1:ll],Modes[1:ll,3],col="red",lwd=2,type="b")
plot(xn$realtime[1:ll],Dist,type="b",ylab="Dist to proj" ,lwd=2,cex.axis=2,cex.lab=2,xlim=c(0,H),ylim=c(0,0.5))
abline(h=0.35,col="grey")
dev.off()

IndexGamma




















