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
xlim=c(0,1)

distance="Lm"
grille=642
		load(paste("RData/XGrid.RData",sep=""))
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )
Gdist<-apply(Gamma,2,GammaDist);

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
	
				NT=NNeighborBlockT(Psitemp)[nOmega] # closest neighbour in Gamma with Lm distance
				Mbelief=Gamma[,NT]
				Mprojprobs<-c(sum(Mbelief[1:nm0]),sum(Mbelief[(nm0+1):(nm0+nm1)]),sum(Mbelief[(nm0+nm1+1):(nm0+nm1+nm2)])) # mass probability over each mode
				dist=distblock(c(Mbelief,Mprojprobs),Psitemp) # Lm distance between filter and its projection
		
				Dist=c(Dist,dist)
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


thresh=0.6
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


pdf(paste("figures/HistDistLm.pdf",sep=""),width=15,height=10)
hist(Dist,30,main="Distances to projected filters",xlim=xlim)
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")


###############################

grille=919
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )
Gdist<-apply(Gamma,2,GammaDist);

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
	
				NT=NNeighborBlockT(Psitemp)[nOmega] # closest neighbour in Gamma with Lm distance
				Mbelief=Gamma[,NT]
				Mprojprobs<-c(sum(Mbelief[1:nm0]),sum(Mbelief[(nm0+1):(nm0+nm1)]),sum(Mbelief[(nm0+nm1+1):(nm0+nm1+nm2)])) # mass probability over each mode
				dist=distblock(c(Mbelief,Mprojprobs),Psitemp) # Lm distance between filter and its projection
		
				Dist=c(Dist,dist)
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


thresh=0.6
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=xlim)
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")
#dev.off()



###############################

grille=1292
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )
Gdist<-apply(Gamma,2,GammaDist);

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
	
				NT=NNeighborBlockT(Psitemp)[nOmega] # closest neighbour in Gamma with Lm distance
				Mbelief=Gamma[,NT]
				Mprojprobs<-c(sum(Mbelief[1:nm0]),sum(Mbelief[(nm0+1):(nm0+nm1)]),sum(Mbelief[(nm0+nm1+1):(nm0+nm1+nm2)])) # mass probability over each mode
				dist=distblock(c(Mbelief,Mprojprobs),Psitemp) # Lm distance between filter and its projection
		
				Dist=c(Dist,dist)
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


thresh=0.55
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=xlim)
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")
#dev.off()


###############################

grille=1346
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )
Gdist<-apply(Gamma,2,GammaDist);

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
	
				NT=NNeighborBlockT(Psitemp)[nOmega] # closest neighbour in Gamma with Lm distance
				Mbelief=Gamma[,NT]
				Mprojprobs<-c(sum(Mbelief[1:nm0]),sum(Mbelief[(nm0+1):(nm0+nm1)]),sum(Mbelief[(nm0+nm1+1):(nm0+nm1+nm2)])) # mass probability over each mode
				dist=distblock(c(Mbelief,Mprojprobs),Psitemp) # Lm distance between filter and its projection
		
				Dist=c(Dist,dist)
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


thresh=0.5
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=xlim)
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")
#dev.off()



###############################

grille=1714
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )
Gdist<-apply(Gamma,2,GammaDist);

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
	
				NT=NNeighborBlockT(Psitemp)[nOmega] # closest neighbour in Gamma with Lm distance
				Mbelief=Gamma[,NT]
				Mprojprobs<-c(sum(Mbelief[1:nm0]),sum(Mbelief[(nm0+1):(nm0+nm1)]),sum(Mbelief[(nm0+nm1+1):(nm0+nm1+nm2)])) # mass probability over each mode
				dist=distblock(c(Mbelief,Mprojprobs),Psitemp) # Lm distance between filter and its projection
		
				Dist=c(Dist,dist)
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


thresh=0.4
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=xlim)
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")
###############################

grille=1807
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )
Gdist<-apply(Gamma,2,GammaDist);

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
	
				NT=NNeighborBlockT(Psitemp)[nOmega] # closest neighbour in Gamma with Lm distance
				Mbelief=Gamma[,NT]
				Mprojprobs<-c(sum(Mbelief[1:nm0]),sum(Mbelief[(nm0+1):(nm0+nm1)]),sum(Mbelief[(nm0+nm1+1):(nm0+nm1+nm2)])) # mass probability over each mode
				dist=distblock(c(Mbelief,Mprojprobs),Psitemp) # Lm distance between filter and its projection
		
				Dist=c(Dist,dist)
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


hist(Dist,30,main="Distances to projected filters",xlim=xlim)
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")


###############################

grille=2144
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )
Gdist<-apply(Gamma,2,GammaDist);

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
	
				NT=NNeighborBlockT(Psitemp)[nOmega] # closest neighbour in Gamma with Lm distance
				Mbelief=Gamma[,NT]
				Mprojprobs<-c(sum(Mbelief[1:nm0]),sum(Mbelief[(nm0+1):(nm0+nm1)]),sum(Mbelief[(nm0+nm1+1):(nm0+nm1+nm2)])) # mass probability over each mode
				dist=distblock(c(Mbelief,Mprojprobs),Psitemp) # Lm distance between filter and its projection
		
				Dist=c(Dist,dist)
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


thresh=0.4
keep=which(Dist>thresh)
prop=round(100*length(keep)/length(Dist),2)


hist(Dist,30,main="Distances to projected filters",xlim=xlim)
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")

###############################

grille=2114
		load(paste("RData/ThetaGrid_",grille,"-",distance,".RData",sep="") )
Gdist<-apply(Gamma,2,GammaDist);

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
	
				NT=NNeighborBlockT(Psitemp)[nOmega] # closest neighbour in Gamma with Lm distance
				Mbelief=Gamma[,NT]
				Mprojprobs<-c(sum(Mbelief[1:nm0]),sum(Mbelief[(nm0+1):(nm0+nm1)]),sum(Mbelief[(nm0+nm1+1):(nm0+nm1+nm2)])) # mass probability over each mode
				dist=distblock(c(Mbelief,Mprojprobs),Psitemp) # Lm distance between filter and its projection
		
				Dist=c(Dist,dist)
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


hist(Dist,30,main="Distances to projected filters",xlim=xlim)
legend("top",legend=paste('Grid size:', grille),bty="n")
legend("right",legend=paste('Proportion:',prop,'%'),bty="n")
abline(v=thresh,col="dark blue")


dev.off()

