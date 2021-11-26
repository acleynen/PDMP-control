set.seed(11)
wd=("R")
source(paste(wd,"SimulationsFunctions.R",sep=""))
source(paste(wd,"Parameters.R",sep=""))
source(paste(wd,"Costs.R",sep=""))
source(paste(wd,"ProgDyn.R",sep=""))
source(paste(wd,"FunctionEvalStrategies.R",sep=""))

flink="id"
trunc=FALSE


nsim=100
cost="simple"

distance="L2"
grille=642
		load(paste("RData/XGrid.RData",sep=""))
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )

Strat1<-NULL
CVs=1; CVi=1
step="all"
	load(paste("RData/ProgDyn_",grille,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))				
FilterReal=c();Dist=c();

IndexGamma<-c()

	for (sim in 1:nsim)
	{
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
		Psin=as.matrix(c(1,rep(0,nOmega-2)))
		if (step==60) Strat=list(dec=rep("non",N),vis=c(60,rep(0,N-1))) else Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(IndexGamma,1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psin<-cbind(Psin,rep(0,nOmega-1)); xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H
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
				FilterReal=cbind(FilterReal,Psitemp)
				NT=NNeighborT(Psitemp)[nOmega]
				Dist=c(Dist,sqrt(sum((Psitemp-Gamma[,NT])^2)))
				Psin<-cbind(Psin,Psitemp)
				indGamma=NT
				IndexGamma<-c(IndexGamma,indGamma)
				Mbelief<-Gamma[,indGamma]
				indvisit=RealTime/delta
				if (cost=="simple")
				{
					if (RealTime<H)
					{
						decision<-Decsimple[indvisit+1,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}	
				} else
				{
					if (RealTime<H)
					{
						decision<-Decint[indvisit,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}					
				}
				yn=c(yn,ynew)
			}
			n<-n+1
		}
	}


thresh=0.35
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


pdf(paste("figures/HistDistL2.pdf",sep=""),width=15,height=10)
hist(Dist,30,main="Distances to projected filters",xlim=c(0,0.6))
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")
#dev.off()


###############################

grille=977
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )

Strat1<-NULL
CVs=1; CVi=1
step="all"
	load(paste("RData/ProgDyn_",grille,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))				

FilterReal=c();Dist=c();
cost="simple"

IndexGamma<-c()

	for (sim in 1:nsim)
	{
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
		Psin=c(1,rep(0,N-1))
		if (step==60) Strat=list(dec=rep("non",N),vis=c(60,rep(0,N-1))) else Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(IndexGamma,1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psinewproj=ngamma; xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H
			} else
			{	
				vmult=which(Time==Strat$vis[n])
				dmult=which(Treatment==Strat$dec[n])
				dind=(vmult-1)*3+dmult
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc) 
					if (ynew<2) ynew=0
				PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Gamma[,Psin[n]]
				eG<-evalGauss(ynew-xval)
				Num=PTemp*eG
				Den=eG%*%PTemp
				Psitemp=Num/as.vector(Den)
				FilterReal=cbind(FilterReal,Psitemp)
				NT=NNeighborT(Psitemp)[nOmega]
				Dist=c(Dist,sqrt(sum((Psitemp-Gamma[,NT])^2)))
				Psin[n+1]<-NT
				indGamma=NT
				IndexGamma<-c(IndexGamma,indGamma)
				Mbelief<-Gamma[,indGamma]
				indvisit=RealTime/delta
				if (cost=="simple")
				{
					if (RealTime<H)
					{
						decision<-Decsimple[indvisit+1,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}	
				} else
				{
					if (RealTime<H)
					{
						decision<-Decint[indvisit,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}					
				}
				yn=c(yn,ynew)
			}
			n<-n+1
		}
	}


thresh=0.35
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=c(0,0.6))
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")
#dev.off()



###############################

grille=1260
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )

Strat1<-NULL
CVs=1; CVi=1
step="all"
	load(paste("RData/ProgDyn_",grille,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))				

FilterReal=c();Dist=c();
cost="simple"

IndexGamma<-c()

	for (sim in 1:nsim)
	{
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
		Psin=c(1,rep(0,N-1))
		if (step==60) Strat=list(dec=rep("non",N),vis=c(60,rep(0,N-1))) else Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(IndexGamma,1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psinewproj=ngamma; xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H
			} else
			{	
				vmult=which(Time==Strat$vis[n])
				dmult=which(Treatment==Strat$dec[n])
				dind=(vmult-1)*3+dmult
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc) 
					if (ynew<2) ynew=0
				PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Gamma[,Psin[n]]
				eG<-evalGauss(ynew-xval)
				Num=PTemp*eG
				Den=eG%*%PTemp
				Psitemp=Num/as.vector(Den)
				FilterReal=cbind(FilterReal,Psitemp)
				NT=NNeighborT(Psitemp)[nOmega]
				Dist=c(Dist,sqrt(sum((Psitemp-Gamma[,NT])^2)))
				Psin[n+1]<-NT
				indGamma=NT
				IndexGamma<-c(IndexGamma,indGamma)
				Mbelief<-Gamma[,indGamma]
				indvisit=RealTime/delta
				if (cost=="simple")
				{
					if (RealTime<H)
					{
						decision<-Decsimple[indvisit+1,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}	
				} else
				{
					if (RealTime<H)
					{
						decision<-Decint[indvisit,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}					
				}
				yn=c(yn,ynew)
			}
			n<-n+1
		}
	}


thresh=0.35
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=c(0,0.6))
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")
#dev.off()


###############################

grille=1212
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )

Strat1<-NULL
CVs=1; CVi=1
step="all"
	load(paste("RData/ProgDyn_",grille,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))				

FilterReal=c();Dist=c();
cost="simple"

IndexGamma<-c()

	for (sim in 1:nsim)
	{
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
		Psin=c(1,rep(0,N-1))
		if (step==60) Strat=list(dec=rep("non",N),vis=c(60,rep(0,N-1))) else Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(IndexGamma,1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psinewproj=ngamma; xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H
			} else
			{	
				vmult=which(Time==Strat$vis[n])
				dmult=which(Treatment==Strat$dec[n])
				dind=(vmult-1)*3+dmult
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc) 
					if (ynew<2) ynew=0
				PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Gamma[,Psin[n]]
				eG<-evalGauss(ynew-xval)
				Num=PTemp*eG
				Den=eG%*%PTemp
				Psitemp=Num/as.vector(Den)
				FilterReal=cbind(FilterReal,Psitemp)
				NT=NNeighborT(Psitemp)[nOmega]
				Dist=c(Dist,sqrt(sum((Psitemp-Gamma[,NT])^2)))
				Psin[n+1]<-NT
				indGamma=NT
				IndexGamma<-c(IndexGamma,indGamma)
				Mbelief<-Gamma[,indGamma]
				indvisit=RealTime/delta
				if (cost=="simple")
				{
					if (RealTime<H)
					{
						decision<-Decsimple[indvisit+1,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}	
				} else
				{
					if (RealTime<H)
					{
						decision<-Decint[indvisit,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}					
				}
				yn=c(yn,ynew)
			}
			n<-n+1
		}
	}


thresh=0.3
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=c(0,0.6))
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")
#dev.off()



###############################

grille=1576
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )

Strat1<-NULL
CVs=1; CVi=1
step="all"
	load(paste("RData/ProgDyn_",grille,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))				

FilterReal=c();Dist=c();
cost="simple"

IndexGamma<-c()

	for (sim in 1:nsim)
	{
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
		Psin=c(1,rep(0,N-1))
		if (step==60) Strat=list(dec=rep("non",N),vis=c(60,rep(0,N-1))) else Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(IndexGamma,1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psinewproj=ngamma; xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H
			} else
			{	
				vmult=which(Time==Strat$vis[n])
				dmult=which(Treatment==Strat$dec[n])
				dind=(vmult-1)*3+dmult
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc) 
					if (ynew<2) ynew=0
				PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Gamma[,Psin[n]]
				eG<-evalGauss(ynew-xval)
				Num=PTemp*eG
				Den=eG%*%PTemp
				Psitemp=Num/as.vector(Den)
				FilterReal=cbind(FilterReal,Psitemp)
				NT=NNeighborT(Psitemp)[nOmega]
				Dist=c(Dist,sqrt(sum((Psitemp-Gamma[,NT])^2)))
				Psin[n+1]<-NT
				indGamma=NT
				IndexGamma<-c(IndexGamma,indGamma)
				Mbelief<-Gamma[,indGamma]
				indvisit=RealTime/delta
				if (cost=="simple")
				{
					if (RealTime<H)
					{
						decision<-Decsimple[indvisit+1,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}	
				} else
				{
					if (RealTime<H)
					{
						decision<-Decint[indvisit,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}					
				}
				yn=c(yn,ynew)
			}
			n<-n+1
		}
	}


thresh=0.3
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=c(0,0.6))
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")


###############################

grille=1582
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )

Strat1<-NULL
CVs=1; CVi=1
step="all"
	load(paste("RData/ProgDyn_",grille,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))				

FilterReal=c();Dist=c();
cost="simple"

IndexGamma<-c()

	for (sim in 1:nsim)
	{
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
		Psin=c(1,rep(0,N-1))
		if (step==60) Strat=list(dec=rep("non",N),vis=c(60,rep(0,N-1))) else Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(IndexGamma,1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psinewproj=ngamma; xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H
			} else
			{	
				vmult=which(Time==Strat$vis[n])
				dmult=which(Treatment==Strat$dec[n])
				dind=(vmult-1)*3+dmult
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc) 
					if (ynew<2) ynew=0
				PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Gamma[,Psin[n]]
				eG<-evalGauss(ynew-xval)
				Num=PTemp*eG
				Den=eG%*%PTemp
				Psitemp=Num/as.vector(Den)
				FilterReal=cbind(FilterReal,Psitemp)
				NT=NNeighborT(Psitemp)[nOmega]
				Dist=c(Dist,sqrt(sum((Psitemp-Gamma[,NT])^2)))
				Psin[n+1]<-NT
				indGamma=NT
				IndexGamma<-c(IndexGamma,indGamma)
				Mbelief<-Gamma[,indGamma]
				indvisit=RealTime/delta
				if (cost=="simple")
				{
					if (RealTime<H)
					{
						decision<-Decsimple[indvisit+1,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}	
				} else
				{
					if (RealTime<H)
					{
						decision<-Decint[indvisit,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}					
				}
				yn=c(yn,ynew)
			}
			n<-n+1
		}
	}


thresh=0.25
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=c(0,0.6))
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")



###############################

grille=2100
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )

Strat1<-NULL
CVs=1; CVi=1
step="all"
	load(paste("RData/ProgDyn_",grille,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))				

FilterReal=c();Dist=c();
cost="simple"

IndexGamma<-c()

	for (sim in 1:nsim)
	{
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
		Psin=c(1,rep(0,N-1))
		if (step==60) Strat=list(dec=rep("non",N),vis=c(60,rep(0,N-1))) else Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(IndexGamma,1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psinewproj=ngamma; xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H
			} else
			{	
				vmult=which(Time==Strat$vis[n])
				dmult=which(Treatment==Strat$dec[n])
				dind=(vmult-1)*3+dmult
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc) 
					if (ynew<2) ynew=0
				PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Gamma[,Psin[n]]
				eG<-evalGauss(ynew-xval)
				Num=PTemp*eG
				Den=eG%*%PTemp
				Psitemp=Num/as.vector(Den)
				FilterReal=cbind(FilterReal,Psitemp)
				NT=NNeighborT(Psitemp)[nOmega]
				Dist=c(Dist,sqrt(sum((Psitemp-Gamma[,NT])^2)))
				Psin[n+1]<-NT
				indGamma=NT
				IndexGamma<-c(IndexGamma,indGamma)
				Mbelief<-Gamma[,indGamma]
				indvisit=RealTime/delta
				if (cost=="simple")
				{
					if (RealTime<H)
					{
						decision<-Decsimple[indvisit+1,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}	
				} else
				{
					if (RealTime<H)
					{
						decision<-Decint[indvisit,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}					
				}
				yn=c(yn,ynew)
			}
			n<-n+1
		}
	}


thresh=0.25
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=c(0,0.6))
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")

###############################

grille=2071
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )

Strat1<-NULL
CVs=1; CVi=1
step="all"
	load(paste("RData/ProgDyn_",grille,"-",distance,"_time",step,"_CV", CVs,".RData",sep=""))				

FilterReal=c();Dist=c();
cost="simple"

IndexGamma<-c()

	for (sim in 1:nsim)
	{
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
		Psin=c(1,rep(0,N-1))
		if (step==60) Strat=list(dec=rep("non",N),vis=c(60,rep(0,N-1))) else Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(IndexGamma,1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psinewproj=ngamma; xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H
			} else
			{	
				vmult=which(Time==Strat$vis[n])
				dmult=which(Treatment==Strat$dec[n])
				dind=(vmult-1)*3+dmult
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc) 
					if (ynew<2) ynew=0
				PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Gamma[,Psin[n]]
				eG<-evalGauss(ynew-xval)
				Num=PTemp*eG
				Den=eG%*%PTemp
				Psitemp=Num/as.vector(Den)
				FilterReal=cbind(FilterReal,Psitemp)
				NT=NNeighborT(Psitemp)[nOmega]
				Dist=c(Dist,sqrt(sum((Psitemp-Gamma[,NT])^2)))
				Psin[n+1]<-NT
				indGamma=NT
				IndexGamma<-c(IndexGamma,indGamma)
				Mbelief<-Gamma[,indGamma]
				indvisit=RealTime/delta
				if (cost=="simple")
				{
					if (RealTime<H)
					{
						decision<-Decsimple[indvisit+1,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}	
				} else
				{
					if (RealTime<H)
					{
						decision<-Decint[indvisit,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}					
				}
				yn=c(yn,ynew)
			}
			n<-n+1
		}
	}


thresh=0.2
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=c(0,0.6))
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")

dev.off()

